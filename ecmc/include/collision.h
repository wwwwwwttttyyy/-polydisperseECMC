#pragma once
#include "cell.h"
#include "geometry.h"
#include <limits>

/**
 * @brief ECMC Collision Event Structure
 * 
 * Describes the nearest collision event when active particle moves along specified direction
 */
struct CollisionEvent {
    double distance;       ///< Event distance (how far to move before collision)
    int targetCellID;      ///< Target particle's cell ID
    int targetSlotK;       ///< Target particle's slot in cell
    
    // Initial relative positions (for subsequent collision geometry calculation)
    double dxInitial;      ///< Initial X direction relative distance
    double dyInitial;      ///< Initial Y direction relative distance
    
    // Pressure contribution (contact distance component at collision moment)
    double deltaIJ;        ///< X direction: sqrt(sigma^2 - dy^2), Y direction: sqrt(sigma^2 - dx^2)
    
    /**
     * @brief Check if valid collision found
     * @return true if collision found within maxDistance range
     */
    bool hasCollision() const {
        return distance < std::numeric_limits<double>::infinity();
    }
    
    /**
     * @brief Get relative X distance at collision moment (for pressure calculation etc.)
     * @note Only valid for X direction movement
     */
    double getDxAtCollision() const {
        return dxInitial - distance;
    }
    
    /**
     * @brief Get relative Y distance at collision moment (for pressure calculation etc.)
     * @note Only valid for Y direction movement
     */
    double getDyAtCollision() const {
        return dyInitial - distance;
    }
    
    /**
     * @brief Default constructor (no collision)
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
 * @brief Search for nearest collision along +X direction
 * 
 * Physical meaning:
 *  - Active particle moves from current position along +X direction
 *  - Other particles remain stationary
 *  - Returns how far to move before first collision occurs
 * 
 * Search strategy:
 *  - Only search forward 6 neighbor cells (directions {2, 5, 8, 1, 4, 7})
 *  - Direction 4 (own cell) requires epsilon check to avoid self-collision
 * 
 * @param grid Spatial indexing structure
 * @param activeCellID Active particle's cell ID
 * @param activeSlotK Active particle's slot in cell
 * @param geom Geometry model (for event distance calculation)
 * @param maxDistance Maximum search distance (usually limited by cell size)
 * @return Nearest collision event (distance = INF if no collision)
 */
CollisionEvent findCollisionX(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance);

/**
 * @brief Search for nearest collision along +Y direction
 * 
 * Physical meaning: Same as findCollisionX, but along +Y direction
 * 
 * Search strategy:
 *  - Only search forward 6 neighbor cells (directions {6, 7, 8, 3, 4, 5})
 * 
 * @param grid Spatial indexing structure
 * @param activeCellID Active particle's cell ID
 * @param activeSlotK Active particle's slot in cell
 * @param geom Geometry model
 * @param maxDistance Maximum search distance
 * @return Nearest collision event
 */
CollisionEvent findCollisionY(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance);
