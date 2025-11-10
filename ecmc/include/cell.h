#pragma once
#include "para.h"
#include "system.h"
#include <vector>
#include <array>
#include <random>

class CellGrid {
public:
    // ========== Construction & Initialization ==========
    /**
     * @brief Construct CellGrid with auto-optimized grid
     * @param sys Particle system
     * @param maxOccupancy Max particles per cell (-1 for auto)
     * 
     * Auto optimization:
     *  - cellSize >= 2.5*radiusMax (safety)
     *  - totalCells ~ sqrt(N) to 2*N (performance balance)
     * 
     * Auto-calculated maxOccupancy (when <= 0):
     *  - Formula: (cellSizeX/(2*rmin) + 1) * (cellSizeY/(2*rmin) + 1) / (sqrt(3)/2) * 1.5
     *  - Meaning: hexagonal close-packing + 50% safety margin
     */
    explicit CellGrid(const System& sys, int maxOccupancy = -1);
    
    /**
     * @brief Construct CellGrid with manual grid specification
     * @param sys Particle system
     * @param nx Number of cells in X
     * @param ny Number of cells in Y
     * @param maxOccupancy Max particles per cell (-1 for auto)
     * 
     * For debugging or special cases, ensure cellSize >= 2*radiusMax
     */
    CellGrid(const System& sys, int nx, int ny, int maxOccupancy = -1);
    
    // ========== Core Operations ==========
    /**
     * @brief Rebuild grid from System
     * @param sys Particle system
     * 
     * Use: initialization or after major position/radius changes
     * Complexity: O(N)
     */
    void rebuild(const System& sys);
    
    /**
     * @brief Randomly choose active particle
     * @param rng Random number generator
     * @return {cellID, slotK} Active particle location
     * 
     * Algorithm: rejection sampling (consistent with HistoricDisks)
     */
    std::pair<int, int> chooseActive(std::mt19937& rng);
    
    /**
     * @brief Update particle position (handles cross-cell moves and syncs global coords)
     * @param sys Particle system (for syncing global coords)
     * @param cellID Current cell ID
     * @param slotK Current slot
     * @param deltaX X displacement
     * @param deltaY Y displacement
     * @return {newCellID, newSlotK} Updated location
     * 
     * Automatically handles:
     *  - Compute new global coords with PBC
     *  - Sync sys.x or sys.y
     *  - Within-cell move: update in-place
     *  - Cross-cell move: swap-tail removal + add to new cell
     */
    std::pair<int, int> moveParticle(System& sys, int cellID, int slotK, 
                                      double deltaX, double deltaY);
    
    // ========== Data Access (read-only, for collision detection) ==========
    double getRelX(int cellID, int slotK) const;
    double getRelY(int cellID, int slotK) const;
    double getRadius(int cellID, int slotK) const;
    int getParticleID(int cellID, int slotK) const;
    int getOccupancy(int cellID) const;
    int getNeighbor(int cellID, int direction) const;
    
    // ========== Generic Neighbor Traversal (Core Interface) ==========
    /**
     * @brief Iterate over neighbor particles in specified directions
     * @param activeCellID Active particle's cell
     * @param activeSlotK Active particle's slot
     * @param directions Direction list to search (nine-grid encoding 0-8)
     * @param callback Callback function with signature:
     *        void(int targetCellID, int targetSlotK, 
     *             double dx, double dy, double targetRadius)
     *        where dx, dy are target's position relative to active (PBC-aware)
     * 
     * Nine-grid encoding:
     *   6  7  8
     *   3  4  5
     *   0  1  2
     * 
     * Usage examples:
     *  - ECMC collision search: directions = {2, 5, 8, 1, 4, 7} (forward 6)
     *  - RDF calculation: directions = {0,1,2,3,4,5,6,7,8} (all neighbors)
     */
    template<typename Func>
    void forEachNeighbor(
        int activeCellID,
        int activeSlotK,
        const std::vector<int>& directions,
        Func callback) const;
    
    /**
     * @brief Get global center coordinates of a cell
     * @param cellID Cell index
     * @return {centerX, centerY} Global coords (box centered at origin)
     */
    std::pair<double, double> getCellCenter(int cellID) const;
    
    // ========== Grid Metadata (read-only) ==========
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    int getTotalCells() const { return totalCells_; }
    double getCellSizeX() const { return cellSizeX_; }
    double getCellSizeY() const { return cellSizeY_; }
    int getMaxOccupancy() const { return maxOccupancy_; }
    
private:
    // ========== Grid Geometry ==========
    int nx_, ny_;                      ///< Grid dimensions
    int totalCells_;                   ///< nx * ny
    double cellSizeX_, cellSizeY_;     ///< Cell dimensions
    int maxOccupancy_;                 ///< Max particles per cell
    std::array<double, 2> boxSize_;    ///< Box size (cached)
    
    // ========== Cell Data (SoA) ==========
    std::vector<int> occupancy_;       ///< Cell occupancy [totalCells]
    std::vector<double> relPosX_;      ///< Relative X coords (flattened)
    std::vector<double> relPosY_;      ///< Relative Y coords (flattened)
    std::vector<double> cellRadius_;   ///< Radii (polydispersity support)
    std::vector<int> particleID_;      ///< Global particle IDs (optional, for debug)
    
    // ========== Neighbor Table (precomputed) ==========
    /// Nine-grid neighbors [totalCells][9]
    /// Encoding: 0-8 correspond to:
    ///   6  7  8
    ///   3  4  5
    ///   0  1  2
    std::vector<std::array<int, 9>> neighbors_;
    
    // ========== Helper Methods ==========
    int cellIndex(int ix, int iy) const;
    std::pair<int, int> cellCoords(int cellID) const;
    void buildNeighborTable();
    void swapRemove(int cellID, int slotK);
    static std::pair<int, int> computeOptimalGrid(
        const std::array<double, 2>& boxSize, 
        double radiusMax, 
        int numParticles);
    
    int flatIndex(int cellID, int slotK) const {
        return cellID * maxOccupancy_ + slotK;
    }
};

// ========== Template Implementation ==========

template<typename Func>
void CellGrid::forEachNeighbor(
    int activeCellID,
    int activeSlotK,
    const std::vector<int>& directions,
    Func callback) const
{
    int activeSlot = flatIndex(activeCellID, activeSlotK);
    double activeRelX = relPosX_[activeSlot];
    double activeRelY = relPosY_[activeSlot];
    
    std::pair<int, int> activeCoords = cellCoords(activeCellID);
    int activeIx = activeCoords.first;
    int activeIy = activeCoords.second;
    double activeCenterX = (activeIx + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
    double activeCenterY = (activeIy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
    
    double activeGlobalX = activeCenterX + activeRelX;
    double activeGlobalY = activeCenterY + activeRelY;
    
    for (int dir : directions) {
        int neighborCellID = neighbors_[activeCellID][dir];
        int ocp = occupancy_[neighborCellID];
        
        std::pair<int, int> neighborCoords = cellCoords(neighborCellID);
        int neighborIx = neighborCoords.first;
        int neighborIy = neighborCoords.second;
        double neighborCenterX = (neighborIx + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
        double neighborCenterY = (neighborIy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
        
        for (int k = 0; k < ocp; ++k) {
            int targetSlot = flatIndex(neighborCellID, k);
            double targetRelX = relPosX_[targetSlot];
            double targetRelY = relPosY_[targetSlot];
            double targetRadius = cellRadius_[targetSlot];
            
            double targetGlobalX = neighborCenterX + targetRelX;
            double targetGlobalY = neighborCenterY + targetRelY;
            
            double dx = targetGlobalX - activeGlobalX;
            double dy = targetGlobalY - activeGlobalY;
            
            dx -= boxSize_[0] * std::round(dx / boxSize_[0]);
            dy -= boxSize_[1] * std::round(dy / boxSize_[1]);
            
            callback(neighborCellID, k, dx, dy, targetRadius);
        }
    }
}
