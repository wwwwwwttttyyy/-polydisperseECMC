#pragma once
#include "cell.h"
#include "geometry.h"
#include <limits>

/**
 * @brief ECMC collision event structure
 * 
 * Describes the nearest collision event when active particle moves along specified direction
 */
struct CollisionEvent {
    double distance;       ///< Event distance (how far to move before collision)
    int targetCellID;      ///< Collision target particle cell ID
    int targetSlotK;       ///< Collision target particle slot in cell
    
    // Initial relative positions (for subsequent collision geometry calculations)
    double dxInitial;      ///< Initial relative distance in x direction
    double dyInitial;      ///< Initial relative distance in y direction
    
    // Pressure contribution (contact distance component at collision)
    double deltaIJ;        ///< X direction: sqrt(σ² - dy²), Y direction: sqrt(σ² - dx²)
    
    /**
     * @brief Check if valid collision found
     * @return true if collision found within maxDistance range
     */
    bool hasCollision() const {
        return distance < std::numeric_limits<double>::infinity();
    }
    
    /**
     * @brief Get relative distance in X direction at collision (for pressure calculation, etc.)
     * @note Only valid for X direction movement
     */
    double getDxAtCollision() const {
        return dxInitial - distance;
    }
    
    /**
     * @brief Get relative distance in Y direction at collision (for pressure calculation, etc.)
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
 * @param grid spatial index structure
 * @param activeCellID active particle cell ID
 * @param activeSlotK active particle slot in cell
 * @param geom geometry model (for calculating event distance)
 * @param maxDistance maximum search distance (usually limited by cell size)
 * @return nearest collision event (if no collision, distance = INF)
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
 * Physical meaning: same as findCollisionX, but along +Y direction
 * 
 * Search strategy:
 *  - Only search forward 6 neighbor cells (directions {6, 7, 8, 3, 4, 5})
 * 
 * @param grid spatial index structure
 * @param activeCellID active particle cell ID
 * @param activeSlotK active particle slot in cell
 * @param geom geometry model
 * @param maxDistance maximum search distance
 * @return nearest collision event
 */
CollisionEvent findCollisionY(
    const CellGrid& grid,
    int activeCellID,
    int activeSlotK,
    const IGeometry& geom,
    double maxDistance);
