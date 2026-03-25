#include "collision.h"
#include <cmath>

// 数值容差（考虑浮点累积误差）
static const double EPSILON = 1e-16;  // 比机器精度大几个数量级

CollisionEvent findCollisionX(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance)
{
    // X 方向搜索：只搜索前方 6 个邻居 cell
    // 九宫格编码：
    //   6  7  8
    //   3  4  5
    //   0  1  2
    // 前方（+X 方向）：2(右下), 5(右中), 8(右上), 1(下), 4(中心), 7(上)
    static const std::vector<int> xDirections = {2, 5, 8, 1, 4, 7};
    
    // 获取活跃粒子信息
    double activeRadius = grid.getRadius(activeCellID, activeSlotK);
    int activeParticleID = grid.getParticleID(activeCellID, activeSlotK);
    
    // 初始化最近碰撞（distance = INF 表示无碰撞）
    CollisionEvent minEvent;
    // minEvent.distance 默认为 INF（由构造函数设置）


    grid.forEachNeighbor(activeCellID, activeSlotK, xDirections,
[&](int targetCellID, int targetSlotK, double dx, double dy, double targetRadius) {

    int targetParticleID = grid.getParticleID(targetCellID, targetSlotK);
    if (activeParticleID == targetParticleID) return;

    if (dx < 1e-12) return;
    double eventDist = geom.XeventDistance(dx, dy, activeRadius, targetRadius);
    if (std::isnan(eventDist)) {
        eventDist = std::numeric_limits<double>::infinity();
    }
    if (eventDist < 0) eventDist = 0.0;
    if (eventDist > maxDistance) return;

    if (eventDist < minEvent.distance) {
        minEvent.distance = eventDist;
        minEvent.targetCellID = targetCellID;
        minEvent.targetSlotK = targetSlotK;
        minEvent.dxInitial = dx;
        minEvent.dyInitial = dy;
        minEvent.deltaIJ = dx - eventDist;
    }
});

    
    return minEvent;
}

CollisionEvent findCollisionY(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance)
{
    // Y 方向搜索：只搜索前方 6 个邻居 cell
    // 前方（+Y 方向）：6(左上), 7(上), 8(右上), 3(左), 4(中心), 5(右)
    static const std::vector<int> yDirections = {6, 7, 8, 3, 4, 5};
    
    // 获取活跃粒子信息
    double activeRadius = grid.getRadius(activeCellID, activeSlotK);
    int activeParticleID = grid.getParticleID(activeCellID, activeSlotK);
    
    // 初始化最近碰撞（distance = INF 表示无碰撞）
    CollisionEvent minEvent;
    
    // 遍历前方邻居
    grid.forEachNeighbor(activeCellID, activeSlotK, yDirections,
        [&](int targetCellID, int targetSlotK, double dx, double dy, double targetRadius) {
            // 避免自碰撞
            int targetParticleID = grid.getParticleID(targetCellID, targetSlotK);
            if (activeParticleID == targetParticleID) {
                return;
            }
            if (dy<1e-12){
                return;
            }
            // 计算事件距离
            double eventDist = geom.YeventDistance(dx, dy, activeRadius, targetRadius);
                if (std::isnan(eventDist)) {
                    eventDist = std::numeric_limits<double>::infinity();
                }
            if(eventDist<0){
                eventDist=0.0;
            }
            if (eventDist > maxDistance) {
                return;
            }

            // 更新最小事件
            if (eventDist < minEvent.distance) {
                minEvent.distance = eventDist;
                minEvent.targetCellID = targetCellID;
                minEvent.targetSlotK = targetSlotK;
                minEvent.dxInitial = dx;
                minEvent.dyInitial = dy;
                minEvent.deltaIJ = dy - eventDist;
            }
        });
    
    return minEvent;
}
