#include "cell.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

// ========== 构造函数 ==========
CellGrid::CellGrid(const System& sys, int maxOccupancy)
    : maxOccupancy_(maxOccupancy),
      boxSize_(sys.boxsize)  // 直接访问 public 成员
{
    // 自动计算最优网格尺寸
    std::pair<int, int> gridSize = computeOptimalGrid(
        boxSize_, 
        sys.radiusMax,  // 直接访问 public 成员
        sys.size()
    );
    
    nx_ = gridSize.first;
    ny_ = gridSize.second;
    totalCells_ = nx_ * ny_;
    cellSizeX_ = boxSize_[0] / nx_;
    cellSizeY_ = boxSize_[1] / ny_;
    
    // 安全性检查
    double minCellSize = std::min(cellSizeX_, cellSizeY_);
    if (minCellSize < 2.0 * sys.radiusMax) {  // 直接访问
        throw std::runtime_error(
            "CellGrid: cellSize < 2*radiusMax! System too dense."
        );
    }
    
    // 初始化数据结构
    occupancy_.resize(totalCells_, 0);
    
   // 计算 maxOccupancy（
    if (maxOccupancy_ <= 0) {
        double rmin = sys.radiusMin > 0 ? sys.radiusMin : sys.radiusMean;
        double sqrtThreeHalf = 0.8660254037844387;  // sqrt(3)/2
        maxOccupancy_ = 1+static_cast<int>(
            (cellSizeX_ / (2.0 * rmin) + 1.0) * 
            (cellSizeY_ / (2.0 * rmin) + 1.0) / 
            sqrtThreeHalf
        );
    }
    
    size_t totalSlots = totalCells_ * maxOccupancy_;
    relPosX_.resize(totalSlots, 0.0);
    relPosY_.resize(totalSlots, 0.0);
    cellRadius_.resize(totalSlots, 0.0);
    particleID_.resize(totalSlots, -1);
    
    neighbors_.resize(totalCells_);
    buildNeighborTable();
}

CellGrid::CellGrid(const System& sys, int nx, int ny, int maxOccupancy)
    : nx_(nx), ny_(ny),
      maxOccupancy_(maxOccupancy),
      boxSize_(sys.boxsize)  // 直接访问 public 成员
{
    totalCells_ = nx_ * ny_;
    cellSizeX_ = boxSize_[0] / nx_;
    cellSizeY_ = boxSize_[1] / ny_;
    
    // 安全性检查: cellSize 必须 > 2 * radiusMax（保证邻居搜索完备）
    double minCellSize = std::min(cellSizeX_, cellSizeY_);
    if (minCellSize < 2.0 * sys.radiusMax) {  // 直接访问
        throw std::runtime_error(
            "CellGrid: cellSize < 2*radiusMax! "
            "Increase nx/ny or decrease system density."
        );
    }
    
    // 初始化数据结构
    occupancy_.resize(totalCells_, 0);
    
    // 智能计算 maxOccupancy（同构造函数 1）
    if (maxOccupancy_ <= 0) {
        double rmin = sys.radiusMin > 0 ? sys.radiusMin : sys.radiusMean;
        double sqrtThreeHalf = 0.8660254037844387;
        maxOccupancy_ = 1 + static_cast<int>(
            (cellSizeX_ / (2.0 * rmin) + 1.0) * 
            (cellSizeY_ / (2.0 * rmin) + 1.0) / 
            sqrtThreeHalf
        );
    }
    
    size_t totalSlots = totalCells_ * maxOccupancy_;
    relPosX_.resize(totalSlots, 0.0);
    relPosY_.resize(totalSlots, 0.0);
    cellRadius_.resize(totalSlots, 0.0);
    particleID_.resize(totalSlots, -1);
    
    neighbors_.resize(totalCells_);
    buildNeighborTable();
}

// ========== 重建网格 ==========
void CellGrid::rebuild(const System& sys) {
    // 清空所有 cell
    std::fill(occupancy_.begin(), occupancy_.end(), 0);
    
    // 遍历所有粒子，插入对应 cell
    for (size_t i = 0; i < sys.size(); ++i) {
        double x = sys.x[i];      // 直接访问 public 成员
        double y = sys.y[i];      // 直接访问 public 成员
        double r = sys.radius[i]; // 直接访问 public 成员
        
        // 全局坐标转 cell 索引（注意：box 中心在原点）
        int ix = static_cast<int>(std::floor((x + boxSize_[0] / 2.0) / cellSizeX_));
        int iy = static_cast<int>(std::floor((y + boxSize_[1] / 2.0) / cellSizeY_));
        
        // PBC wrap
        ix = (ix + nx_) % nx_;
        iy = (iy + ny_) % ny_;
        
        int cellID = cellIndex(ix, iy);
        int ocp = occupancy_[cellID];
        
        if (ocp >= maxOccupancy_) {
            throw std::runtime_error(
                "CellGrid::rebuild: Cell overflow! "
                "Increase maxOccupancy or adjust grid."
            );
        }
        
        // Cell 中心全局坐标
        double centerX = (ix + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
        double centerY = (iy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
        
        // 计算相对坐标（考虑 PBC）
        double relX = x - centerX;
        double relY = y - centerY;
        
        // PBC 校正（确保相对坐标在 [-cellSize/2, cellSize/2]）
        relX -= std::round(relX / cellSizeX_) * cellSizeX_;
        relY -= std::round(relY / cellSizeY_) * cellSizeY_;
        
        // 存储
        int slot = flatIndex(cellID, ocp);
        relPosX_[slot] = relX;
        relPosY_[slot] = relY;
        cellRadius_[slot] = r;
        particleID_[slot] = static_cast<int>(i);
        
        occupancy_[cellID]++;
    }
}

// ========== 选择活跃粒子 ==========
std::pair<int, int> CellGrid::chooseActive(std::mt19937& rng) {
    std::uniform_int_distribution<int> cellDist(0, totalCells_ - 1);
    
    while (true) {
        int cellID = cellDist(rng);
        int ocp = occupancy_[cellID];
        
        if (ocp > 0) {
            std::uniform_int_distribution<int> slotDist(0, ocp - 1);
            int slotK = slotDist(rng);
            return std::make_pair(cellID, slotK);
        }
    }
}

// ========== 移动粒子 ==========
std::pair<int, int> CellGrid::moveParticle(
    System& sys, int cellID, int slotK, double deltaX, double deltaY)
{
    int slot = flatIndex(cellID, slotK);
    double newRelX = relPosX_[slot] + deltaX;
    double newRelY = relPosY_[slot] + deltaY;
    
    // 获取粒子 ID
    int particleID = particleID_[slot];
    
    // 检查是否跨 cell
    bool crossX = std::abs(newRelX) > cellSizeX_ / 2.0;
    bool crossY = std::abs(newRelY) > cellSizeY_ / 2.0;
    
    if (!crossX && !crossY) {
        // 在 cell 内移动，直接更新
        relPosX_[slot] = newRelX;
        relPosY_[slot] = newRelY;
        
        // 同步更新 System 全局坐标
        if (deltaX != 0.0) {
            sys.x[particleID] += deltaX;
            sys.x[particleID] -= boxSize_[0] * std::round(sys.x[particleID] / boxSize_[0]);
        }
        if (deltaY != 0.0) {
            sys.y[particleID] += deltaY;
            sys.y[particleID] -= boxSize_[1] * std::round(sys.y[particleID] / boxSize_[1]);
        }
        
        return std::make_pair(cellID, slotK);
    }
    
    // 跨 cell：需要重新分配
    std::pair<int, int> coords = cellCoords(cellID);
    int ix = coords.first;
    int iy = coords.second;
    
    // 根据相对坐标判断跨越方向 (利用九宫格邻居表)
    int dx_dir = 0;
    int dy_dir = 0;
    
    if (newRelX > cellSizeX_ / 2.0) dx_dir = 1;
    else if (newRelX < -cellSizeX_ / 2.0) dx_dir = -1;
    
    if (newRelY > cellSizeY_ / 2.0) dy_dir = 1;
    else if (newRelY < -cellSizeY_ / 2.0) dy_dir = -1;
    
    // 九宫格编码: dir = (dy + 1) * 3 + (dx + 1)
    int neighborDir = (dy_dir + 1) * 3 + (dx_dir + 1);
    int newCellID = neighbors_[cellID][neighborDir];
    
    // 新 cell 中心
    std::pair<int, int> newCoords = cellCoords(newCellID);
    int newIx = newCoords.first;
    int newIy = newCoords.second;
    double newCenterX = (newIx + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
    double newCenterY = (newIy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
    
    // 原 cell 中心
    double oldCenterX = (ix + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
    double oldCenterY = (iy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
    
    // 全局坐标(用于计算新相对坐标)
    double globalX = oldCenterX + relPosX_[slot] + deltaX;
    double globalY = oldCenterY + relPosY_[slot] + deltaY;
    
    // 应用 PBC 到全局坐标
    globalX -= boxSize_[0] * std::round(globalX / boxSize_[0]);
    globalY -= boxSize_[1] * std::round(globalY / boxSize_[1]);
    
    // 同步更新 System 全局坐标
    sys.x[particleID] = globalX;
    sys.y[particleID] = globalY;
    
    // 新相对坐标
    double finalRelX = globalX - newCenterX;
    double finalRelY = globalY - newCenterY;
    
    // PBC 校正
    finalRelX -= std::round(finalRelX / cellSizeX_) * cellSizeX_;
    finalRelY -= std::round(finalRelY / cellSizeY_) * cellSizeY_;
    
    // 保存粒子属性
    double radius = cellRadius_[slot];
    int particleIdx = particleID_[slot];
    
    // Swap-tail 删除（O(1)）
    swapRemove(cellID, slotK);
    
    // 添加到新 cell
    int newOcp = occupancy_[newCellID];
    if (newOcp >= maxOccupancy_) {
        throw std::runtime_error("CellGrid::moveParticle: New cell overflow!");
    }
    
    int newSlot = flatIndex(newCellID, newOcp);
    relPosX_[newSlot] = finalRelX;
    relPosY_[newSlot] = finalRelY;
    cellRadius_[newSlot] = radius;
    particleID_[newSlot] = particleIdx;
    
    occupancy_[newCellID]++;
    
    return std::make_pair(newCellID, newOcp);
}

// ========== 数据访问 ==========
double CellGrid::getRelX(int cellID, int slotK) const {
    return relPosX_[flatIndex(cellID, slotK)];
}

double CellGrid::getRelY(int cellID, int slotK) const {
    return relPosY_[flatIndex(cellID, slotK)];
}

double CellGrid::getRadius(int cellID, int slotK) const {
    return cellRadius_[flatIndex(cellID, slotK)];
}

int CellGrid::getParticleID(int cellID, int slotK) const {
    return particleID_[flatIndex(cellID, slotK)];
}

int CellGrid::getOccupancy(int cellID) const {
    return occupancy_[cellID];
}

int CellGrid::getNeighbor(int cellID, int direction) const {
    return neighbors_[cellID][direction];
}

std::pair<double, double> CellGrid::getCellCenter(int cellID) const {
    std::pair<int, int> coords = cellCoords(cellID);
    int ix = coords.first;
    int iy = coords.second;
    
    double centerX = (ix + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
    double centerY = (iy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
    
    return std::make_pair(centerX, centerY);
}

// ========== 辅助方法 ==========
int CellGrid::cellIndex(int ix, int iy) const {
    return iy * nx_ + ix;
}

std::pair<int, int> CellGrid::cellCoords(int cellID) const {
    int iy = cellID / nx_;
    int ix = cellID % nx_;
    return std::make_pair(ix, iy);
}

void CellGrid::buildNeighborTable() {
    /**
     * 九宫格编码:
     *   6  7  8
     *   3  4  5
     *   0  1  2
     * 
     * 对应偏移:
     *   (-1,+1) (0,+1) (+1,+1)
     *   (-1, 0) (0, 0) (+1, 0)
     *   (-1,-1) (0,-1) (+1,-1)
     */
    const std::array<std::pair<int, int>, 9> offsets = {{
        {-1, -1}, {0, -1}, {1, -1},
        {-1,  0}, {0,  0}, {1,  0},
        {-1,  1}, {0,  1}, {1,  1}
    }};
    
    for (int iy = 0; iy < ny_; ++iy) {
        for (int ix = 0; ix < nx_; ++ix) {
            int cellID = cellIndex(ix, iy);
            
            for (int dir = 0; dir < 9; ++dir) {
                int dx = offsets[dir].first;
                int dy = offsets[dir].second;
                
                int nix = (ix + dx + nx_) % nx_;  // PBC
                int niy = (iy + dy + ny_) % ny_;
                
                neighbors_[cellID][dir] = cellIndex(nix, niy);
            }
        }
    }
}

void CellGrid::swapRemove(int cellID, int slotK) {
    int ocp = occupancy_[cellID];
    
    if (slotK >= ocp || ocp == 0) {
        throw std::out_of_range("CellGrid::swapRemove: Invalid slot!");
    }
    
    // Swap with last element
    int lastSlot = flatIndex(cellID, ocp - 1);
    int currSlot = flatIndex(cellID, slotK);
    
    if (slotK != ocp - 1) {
        relPosX_[currSlot] = relPosX_[lastSlot];
        relPosY_[currSlot] = relPosY_[lastSlot];
        cellRadius_[currSlot] = cellRadius_[lastSlot];
        particleID_[currSlot] = particleID_[lastSlot];
    }
    
    occupancy_[cellID]--;
}

std::pair<int, int> CellGrid::computeOptimalGrid(
    const std::array<double, 2>& boxSize,
    double radiusMax,
    int numParticles)
{
    /**
     * 自动网格尺寸计算策略:
     *  1. 确保 cellSize >= 2.5 * radiusMax（安全余量）
     *  2. 优化 cell 数量：约 N^0.5 cells（经验公式）
     * 
     * 算法: 从最小 cellSize 开始，调整到合理粒子密度
     */
    double minCellSize = 2.5 * radiusMax;
    
    int nx = std::max(1, static_cast<int>(boxSize[0] / minCellSize));
    int ny = std::max(1, static_cast<int>(boxSize[1] / minCellSize));
    
    // 避免过度细分（经验上限：totalCells ~ 2*N）
    int totalCells = nx * ny;
    int maxCells = std::max(10, 2 * numParticles);
    
    if (totalCells > maxCells) {
        double scale = std::sqrt(static_cast<double>(maxCells) / totalCells);
        nx = std::max(1, static_cast<int>(nx * scale));
        ny = std::max(1, static_cast<int>(ny * scale));
    }
    
    return std::make_pair(nx, ny);
}
