#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>

class System {
public:
    // ========== 构造与初始化 ==========
    System() = default;
    explicit System(size_t n, const std::array<double, 2>& box_);  ///< 预分配 n 个粒子并设置盒边界
    
    std::vector<double> x;        ///< x 坐标 [N]
    std::vector<double> y;        ///< y 坐标 [N]
    std::vector<double> radius;   ///< 半径 [N]（多分散关键！）
    std::vector<int> typeID;      ///< 类型 ID [N]（可选，用于统计）
    std::array<double, 2> boxsize;    ///< {Lx, Ly}，应由外部根据 PhysicalPara 设置，在实际模拟中假设盒子中心在原点
    // 储存最大半径用于构建cell list
    double radiusMax = 0.0;       ///< max(radius)
    double radiusMean = 0.0;      ///< mean(radius)
    double radiusMin = std::numeric_limits<double>::infinity();  ///< min(radius)

    size_t size() const { return x.size(); }
    
    // 预分配内存
    void reserve(size_t n);

    // 清空所有数据，释放空间
    void clear();
    
    // ========== 粒子操作 ==========
    void addParticle(double x_, double y_, double r_, int tid = 0);
    void removeParticle(size_t index);
    void updateMaxMeanRadius();
    
    // ========== 文件 IO ==========
    /**
     * @brief 保存系统构型到文本文件（人类可读）
     * @param filename 文件名
     * 
     * 格式：
     *   第1行: Lx Ly
     *   后续行: x y radius
     */
    void saveConfig(const std::string& filename) const;
    
    /**
     * @brief 从文本文件加载系统构型
     * @param filename 文件名
     * 
     * 格式与 saveConfig 相同
     * 会清空当前系统并加载新数据
     */
    void loadConfig(const std::string& filename);
    
    /**
     * @brief 保存系统构型到二进制文件（高效）
     * @param filename 文件名
     * 
     * 格式（所有数据为 double）：
     *   - Lx, Ly (盒子尺寸)
     *   - N (粒子数，作为 double 存储)
     *   - x[0], y[0], r[0]
     *   - x[1], y[1], r[1]
     *   - ...
     */
    void saveConfigBinary(const std::string& filename) const;
    
    /**
     * @brief 从二进制文件加载系统构型
     * @param filename 文件名
     */
    void loadConfigBinary(const std::string& filename);
    
    // 盒子操作与global PBC实现
    //初始位型重合检查.false:无重合.true:有重合
    bool checkOverlap(double overlapThreshold = 1e-9) const;
    //拉回主盒
    void applyPBC(size_t i);
    //global PBC
    std::array<double, 2> displacement(size_t i, size_t j) const;
    //广义距离
    double distance(size_t i, size_t j) const;
    
    // ========== 统计量计算 ==========
    /**
     * @brief 计算系统的填充分数 (packing fraction)
     * @return φ = Σ(πr²) / (Lx * Ly)
     */
    double getPackingFraction() const;
};
