#pragma once
#include "para.h"
#include "system.h"
#include <vector>
#include <array>
#include <random>
#include "geometry.h"
class CellGrid {
public:
    // ========== 构造与初始化 ==========
    /**
     * @brief 构造 CellGrid（自动计算最优网格）
     * @param sys 粒子系统（用于获取 box、radiusMax、radiusMin、N）
     * @param maxOccupancy 单 cell 最大容量（默认 -1 自动计算）
     * 
     * 网格尺寸自动优化：
     *  - cellSize >= 2.5 * radiusMax（安全边界）
     *  - totalCells ~ sqrt(N) 到 2*N（性能平衡）
     * 
     * maxOccupancy 自动计算（当 <= 0 时）：
     *  - 公式: (cellSizeX/(2*rmin) + 1) * (cellSizeY/(2*rmin) + 1) / (sqrt(3)/2) * 1.5
     *  - 物理意义: 六角密堆积理论值 + 50% 安全余量
     */
    explicit CellGrid(const System& sys, int maxOccupancy = -1);
    
    /**
     * @brief 构造 CellGrid（手动指定网格）
     * @param sys 粒子系统
     * @param nx X 方向 cell 数量
     * @param ny Y 方向 cell 数量
     * @param maxOccupancy 单 cell 最大容量（-1 自动计算，否则使用指定值）
     * 
     * 用于调试或特殊场景，需确保 cellSize >= 2*radiusMax
     */
    CellGrid(const System& sys, int nx, int ny, int maxOccupancy = -1);
    
    // ========== 核心操作 ==========
    /**
     * @brief 从 System 重建整个网格
     * @param sys 粒子系统
     * 
     * 用途: 初始化或粒子位置/半径大幅变化后调用
     * 复杂度: O(N)
     */
    void rebuild(const System& sys);
    
    /**
     * @brief 随机选择活跃粒子
     * @param rng 随机数生成器
     * @return {cellID, slotK} 返回活跃粒子在 cell 中的位置
     * 
     * 算法: 拒绝采样（与 HistoricDisks 一致）
     */
    std::pair<int, int> chooseActive(std::mt19937& rng);
    
    /**
     * @brief 更新粒子位置（处理跨 cell 移动并同步 System 全局坐标）
     * @param sys 粒子系统(用于同步全局坐标)
     * @param cellID 当前 cell ID
     * @param slotK 当前 slot
     * @param deltaX X 方向位移
     * @param deltaY Y 方向位移
     * @return {newCellID, newSlotK} 更新后的位置
     * 
     * 自动处理:
     *  - 计算新全局坐标并应用 PBC
     *  - 同步更新 sys.x 或 sys.y
     *  - 在 cell 内移动: 就地更新相对坐标
     *  - 跨 cell 移动: swap-tail 删除 + 新 cell 末尾添加
     */
    std::pair<int, int> moveParticle(System& sys, int cellID, int slotK, 
                                      double deltaX, double deltaY);
    
    // ========== 数据访问（只读，供碰撞检测使用）==========
    /// 获取相对坐标
    double getRelX(int cellID, int slotK) const;
    double getRelY(int cellID, int slotK) const;
    
    /// 获取半径
    double getRadius(int cellID, int slotK) const;
    
    /// 获取全局粒子 ID（可选，用于调试）
    int getParticleID(int cellID, int slotK) const;
    
    /// 获取 cell 占据数
    int getOccupancy(int cellID) const;
    
    /// 获取邻居 cell ID（九宫格索引 0-8）
    int getNeighbor(int cellID, int direction) const;
    
    // ========== 通用邻居遍历（核心接口）==========
    /**
     * @brief 遍历指定方向的邻居粒子
     * @param activeCellID 活跃粒子所在 cell
     * @param activeSlotK 活跃粒子 slot
     * @param directions 要搜索的方向列表（九宫格编码 0-8）
     * @param callback 回调函数，签名：
     *        void(int targetCellID, int targetSlotK, 
     *             double dx, double dy, double targetRadius)
     *        其中 dx, dy 是目标粒子相对于活跃粒子的位置（已处理 PBC）
     * 
     * 九宫格编码：
     *   6  7  8
     *   3  4  5
     *   0  1  2
     * 
     * 用途示例：
     *  - ECMC 碰撞搜索：directions = {2, 5, 8, 1, 4, 7} (前方 6 个)
     *  - RDF 计算：directions = {0,1,2,3,4,5,6,7,8} (所有邻居)
     */
    template<typename Func>
    void forEachNeighbor(
        int activeCellID,
        int activeSlotK,
        const std::vector<int>& directions,
        Func callback) const;
    
    /**
     * @brief 获取 cell 的全局中心坐标
     * @param cellID cell 索引
     * @return {centerX, centerY} 全局坐标（盒子中心在原点）
     */
    std::pair<double, double> getCellCenter(int cellID) const;
    
    // ========== 网格元数据（只读）==========
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    int getTotalCells() const { return totalCells_; }
    double getCellSizeX() const { return cellSizeX_; }
    double getCellSizeY() const { return cellSizeY_; }
    int getMaxOccupancy() const { return maxOccupancy_; }
    
private:
    // ========== 网格几何参数 ==========
    int nx_, ny_;                      ///< 网格维度
    int totalCells_;                   ///< nx * ny
    double cellSizeX_, cellSizeY_;     ///< Cell 边长
    int maxOccupancy_;                 ///< 单 cell 最大容量
    std::array<double, 2> boxSize_;    ///< 盒边界（缓存）
    const System* sysPtr_;             ///< 指向系统（用于热区坐标访问）
    
    // ========== Cell 数据（SoA in cell）==========
    /// Cell 占据数 [totalCells]
    std::vector<int> occupancy_;
    
    /// 相对坐标 [totalCells][maxOccupancy]（展平存储）
    std::vector<double> relPosX_;
    std::vector<double> relPosY_;
    
    /// 半径 [totalCells][maxOccupancy]（多分散关键！）
    std::vector<double> cellRadius_;
    
    /// 全局粒子 ID [totalCells][maxOccupancy]（可选，调试用）
    std::vector<int> particleID_;
    
    // ========== 邻居表（预计算）==========
    /// 九宫格邻居 [totalCells][9]
    /// 编码: 0-8 对应:
    ///   6  7  8
    ///   3  4  5
    ///   0  1  2
    std::vector<std::array<int, 9>> neighbors_;
    
    // ========== 辅助方法 ==========
    /// Cell ID 映射: (ix, iy) -> cellID
    int cellIndex(int ix, int iy) const;
    
    /// Cell ID 反映射: cellID -> (ix, iy)
    std::pair<int, int> cellCoords(int cellID) const;
    
    /// 构建邻居表（考虑 PBC）
    void buildNeighborTable();
    
    /// Swap-tail 删除（O(1)）
    void swapRemove(int cellID, int slotK);
    
    /// 自动计算最优网格尺寸
    static std::pair<int, int> computeOptimalGrid(
        const std::array<double, 2>& boxSize, 
        double radiusMax, 
        int numParticles);
    
    /// 二维索引展平: cellID * maxOccupancy + slotK
    int flatIndex(int cellID, int slotK) const {
         return static_cast<size_t>(cellID) * maxOccupancy_ + slotK;
    }
};

// ========== Template 方法实现 ==========

template<typename Func>
void CellGrid::forEachNeighbor(
    int activeCellID,
    int activeSlotK,
    const std::vector<int>& directions,
    Func callback) const
{
    // Fast path: use System global coords; skip per-neighbor cell-center + floor/round work
    const int activeSlot = flatIndex(activeCellID, activeSlotK);
    const int activePID = particleID_[activeSlot];
    const double ax = sysPtr_->x[activePID];
    const double ay = sysPtr_->y[activePID];

    const double halfBoxX = 0.5 * boxSize_[0];
    const double halfBoxY = 0.5 * boxSize_[1];
    auto wrap = [](double d, double L, double halfL) {
        if (d > halfL)       d -= L;
        else if (d < -halfL) d += L;
        return d;
    };

    for (int dir : directions) {
        int neighborCellID = neighbors_[activeCellID][dir];
        int ocp = occupancy_[neighborCellID];

        for (int k = 0; k < ocp; ++k) {
            int targetSlot = flatIndex(neighborCellID, k);
            int targetPID = particleID_[targetSlot];

            double dx = sysPtr_->x[targetPID] - ax;
            double dy = sysPtr_->y[targetPID] - ay;
            dx = wrap(dx, boxSize_[0], halfBoxX);
            dy = wrap(dy, boxSize_[1], halfBoxY);

            double targetRadius = cellRadius_[targetSlot];
            callback(neighborCellID, k, dx, dy, targetRadius);
        }
    }
}

