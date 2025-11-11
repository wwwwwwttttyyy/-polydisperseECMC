#include "collision.h"
#include <cmath>

// Numerical tolerance (considering floating-point accumulation error)
static const double EPSILON = 1e-12;  // Several orders of magnitude larger than machine precision

CollisionEvent findCollisionX(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance)
{
    // X direction search: only search forward 6 neighbor cells
    // 3x3 grid encoding:
    //   6  7  8
    //   3  4  5
    //   0  1  2
    // Forward (+X direction): 2(bottom-right), 5(right-middle), 8(top-right), 1(bottom), 4(center), 7(top)
    static const std::vector<int> xDirections = {2, 5, 8, 1, 4, 7};
    
    // Get active particle information
    double activeRadius = grid.getRadius(activeCellID, activeSlotK);
    int activeParticleID = grid.getParticleID(activeCellID, activeSlotK);
    
    // Initialize nearest collision (distance = INF means no collision)
    CollisionEvent minEvent;
    // minEvent.distance defaults to INF (set by constructor)
    
    // Iterate over forward neighbors
    grid.forEachNeighbor(activeCellID, activeSlotK, xDirections,
        [&](int targetCellID, int targetSlotK, double dx, double dy, double targetRadius) {
            // 1. Avoid self-collision (same particle)
            int targetParticleID = grid.getParticleID(targetCellID, targetSlotK);
            if (activeParticleID == targetParticleID) {
                return;
            }
            
            // 2. Quick rejection: dx must be positive (target is forward)
            // EPSILON avoids self-collision due to floating-point errors
            if (dx <= EPSILON) {
                return;
            }
            
            // 3. Calculate event distance
            double eventDist = geom.XeventDistance(dx, dy, activeRadius, targetRadius);
            
            // 4. Robustness check: detect numerical errors
            if (std::isnan(eventDist)) {
                return;  // NaN indicates calculation error, skip
            }
            
            if (eventDist < 0.0) {
                return;  // Negative distance is physical error (particles already overlapping), skip
            }
            
            // 5. INF means never colliding (too far in vertical direction)
            if (std::isinf(eventDist)) {
                return;  // Normal case, skip this particle
            }
            
            // 6. Distance limit
            if (eventDist > maxDistance) {
                return;
            }
            
            // 7. Update minimum event
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
    // Y direction search: only search forward 6 neighbor cells
    // Forward (+Y direction): 6(top-left), 7(top), 8(top-right), 3(left), 4(center), 5(right)
    static const std::vector<int> yDirections = {6, 7, 8, 3, 4, 5};
    
    // Get active particle information
    double activeRadius = grid.getRadius(activeCellID, activeSlotK);
    int activeParticleID = grid.getParticleID(activeCellID, activeSlotK);
    
    // Initialize nearest collision (distance = INF means no collision)
    CollisionEvent minEvent;
    
    // Iterate over forward neighbors
    grid.forEachNeighbor(activeCellID, activeSlotK, yDirections,
        [&](int targetCellID, int targetSlotK, double dx, double dy, double targetRadius) {
            // 1. Avoid self-collision
            int targetParticleID = grid.getParticleID(targetCellID, targetSlotK);
            if (activeParticleID == targetParticleID) {
                return;
            }
            
            // 2. Quick rejection: dy must be positive (target is forward)
            if (dy <= EPSILON) {
                return;
            }
            
            // 3. Calculate event distance
            double eventDist = geom.YeventDistance(dx, dy, activeRadius, targetRadius);
            
            // 4. Robustness check: detect numerical errors
            if (std::isnan(eventDist)) {
                return;  // NaN indicates calculation error
            }
            
            if (eventDist < 0.0) {
                return;  // Negative distance is physical error
            }
            
            // 5. INF means never colliding
            if (std::isinf(eventDist)) {
                return;  // Normal case, skip
            }
            
            // 6. Distance limit
            if (eventDist > maxDistance) {
                return;
            }
            
            // 7. Update minimum event
            if (eventDist < minEvent.distance) {
                minEvent.distance = eventDist;
                minEvent.targetCellID = targetCellID;
                minEvent.targetSlotK = targetSlotK;
                minEvent.dxInitial = dx;
                minEvent.dyInitial = dy;
                // deltaIJ = dy - eventDist (Y-direction contact distance at collision time)
                minEvent.deltaIJ = dy - eventDist;
            }
        });
    
    return minEvent;
}
