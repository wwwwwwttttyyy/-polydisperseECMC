#pragma once
#include "cell.h"
#include "geometry.h"
#include <limits>

/**
 * @brief ECMC 碰撞事件结构
 * 
 * 描述活跃粒子沿指定方向移动时的最近碰撞事件
 */
struct CollisionEvent {
    double distance;       ///< 事件距离（移动多远才碰撞）
    int targetCellID;      ///< 碰撞目标粒子所在 cell ID
    int targetSlotK;       ///< 碰撞目标粒子在 cell 中的 slot
    
    // 初始相对位置（用于后续计算碰撞几何）
    double dxInitial;      ///< 初始时刻 x 方向的相对距离
    double dyInitial;      ///< 初始时刻 y 方向的相对距离
    
    // 压强贡献（碰撞时刻的接触距离分量）
    double deltaIJ;        ///< X方向：sqrt(σ² - dy²)，Y方向：sqrt(σ² - dx²)
    
    /**
     * @brief 判断是否找到有效碰撞
     * @return true 如果在 maxDistance 范围内找到碰撞
     */
    bool hasCollision() const {
        return distance < std::numeric_limits<double>::infinity();
    }
    
    /**
     * @brief 获取碰撞时刻 X 方向的相对距离（用于压强计算等）
     * @note 仅在 X 方向移动时有效
     */
    double getDxAtCollision() const {
        return dxInitial - distance;
    }
    
    /**
     * @brief 获取碰撞时刻 Y 方向的相对距离（用于压强计算等）
     * @note 仅在 Y 方向移动时有效
     */
    double getDyAtCollision() const {
        return dyInitial - distance;
    }
    
    /**
     * @brief 默认构造（无碰撞）
     */
    CollisionEvent()
        : distance(std::numeric_limits<double>::infinity()),
          targetCellID(-1),
          targetSlotK(-1),
          dxInitial(0.0),
          dyInitial(0.0),
          deltaIJ(0.0) {}
};

/**
 * @brief 搜索沿 +X 方向的最近碰撞
 * 
 * 物理含义：
 *  - 活跃粒子从当前位置沿 +X 方向移动
 *  - 其他粒子保持静止
 *  - 返回移动多远后会发生第一次碰撞
 * 
 * 搜索策略：
 *  - 只搜索前方 6 个邻居 cell（方向 {2, 5, 8, 1, 4, 7}）
 *  - 方向 4（自己所在 cell）需要通过 epsilon 检查避免自碰撞
 * 
 * @param grid 空间索引结构
 * @param activeCellID 活跃粒子所在 cell ID
 * @param activeSlotK 活跃粒子在 cell 中的 slot
 * @param geom 几何模型（用于计算事件距离）
 * @param maxDistance 最大搜索距离（通常受 cell 大小限制）
 * @return 最近的碰撞事件（如果无碰撞，distance = INF）
 */
CollisionEvent findCollisionX(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance);

/**
 * @brief 搜索沿 +Y 方向的最近碰撞
 * 
 * 物理含义：同 findCollisionX，但沿 +Y 方向
 * 
 * 搜索策略：
 *  - 只搜索前方 6 个邻居 cell（方向 {6, 7, 8, 3, 4, 5}）
 * 
 * @param grid 空间索引结构
 * @param activeCellID 活跃粒子所在 cell ID
 * @param activeSlotK 活跃粒子在 cell 中的 slot
 * @param geom 几何模型
 * @param maxDistance 最大搜索距离
 * @return 最近的碰撞事件
 */
CollisionEvent findCollisionY(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance);
