#pragma once
#include "para.h"
#include "system.h"
#include <vector>
#include <array>
#include <random>

class CellGrid {
public:
    // ========== Construction and Initialization ==========
    /**
     * @brief Construct CellGrid (automatically computes optimal grid)
     * @param sys particle system (for box, radiusMax, radiusMin, N)
     * @param maxOccupancy maximum capacity per cell (default -1 auto-computed)
     * 
     * Grid size auto-optimization:
     *  - cellSize >= 2.5 * radiusMax (safety margin)
     *  - totalCells ~ sqrt(N) to 2*N (performance balance)
     * 
     * maxOccupancy auto-computed (when <= 0):
     *  - Formula: (cellSizeX/(2*rmin) + 1) * (cellSizeY/(2*rmin) + 1) / (sqrt(3)/2) * 1.5
     *  - Physical meaning: hexagonal close packing theoretical value + 50% safety margin
     */
    explicit CellGrid(const System& sys, int maxOccupancy = -1);
    
    /**
     * @brief Construct CellGrid (manually specified grid)
     * @param sys particle system
     * @param nx number of cells in X direction
     * @param ny number of cells in Y direction
     * @param maxOccupancy maximum capacity per cell (-1 auto-computed, otherwise use specified value)
     * 
     * For debugging or special scenarios, ensure cellSize >= 2*radiusMax
     */
    CellGrid(const System& sys, int nx, int ny, int maxOccupancy = -1);
    
    // ========== Core Operations ==========
    /**
     * @brief Rebuild entire grid from System
     * @param sys particle system
     * 
     * Use: initialization or after significant particle position/radius changes
     * Complexity: O(N)
     */
    void rebuild(const System& sys);
    
    /**
     * @brief Randomly select active particle
     * @param rng random number generator
     * @return {cellID, slotK} active particle position in cell
     * 
     * Algorithm: rejection sampling (consistent with HistoricDisks)
     */
    std::pair<int, int> chooseActive(std::mt19937& rng);
    
    /**
     * @brief Update particle position (handles cross-cell movement and syncs System global coordinates)
     * @param sys particle system (for global coordinate sync)
     * @param cellID current cell ID
     * @param slotK current slot
     * @param deltaX X displacement
     * @param deltaY Y displacement
     * @return {newCellID, newSlotK} updated position
     * 
     * Automatically handles:
     *  - Compute new global coordinates and apply PBC
     *  - Synchronously update sys.x or sys.y
     *  - In-cell movement: update relative coordinates in place
     *  - Cross-cell movement: swap-tail delete + add to new cell end
     */
    std::pair<int, int> moveParticle(System& sys, int cellID, int slotK, 
                                      double deltaX, double deltaY);
    
    // ========== Data Access (read-only, for collision detection) ==========
    /// Get relative coordinates
    double getRelX(int cellID, int slotK) const;
    double getRelY(int cellID, int slotK) const;
    
    /// Get radius
    double getRadius(int cellID, int slotK) const;
    
    /// Get global particle ID (optional, for debugging)
    int getParticleID(int cellID, int slotK) const;
    
    /// Get cell occupancy
    int getOccupancy(int cellID) const;
    
    /// Get neighbor cell ID (3x3 grid index 0-8)
    int getNeighbor(int cellID, int direction) const;
    
    // ========== Generic Neighbor Traversal (core interface) ==========
    /**
     * @brief Iterate over neighbor particles in specified directions
     * @param activeCellID active particle cell
     * @param activeSlotK active particle slot
     * @param directions list of directions to search (3x3 grid encoding 0-8)
     * @param callback callback function, signature:
     *        void(int targetCellID, int targetSlotK, 
     *             double dx, double dy, double targetRadius)
     *        where dx, dy are target particle positions relative to active particle (PBC handled)
     * 
     * 3x3 grid encoding:
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
     * @brief Get cell's global center coordinates
     * @param cellID cell index
     * @return {centerX, centerY} global coordinates (box center at origin)
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
    // ========== Grid Geometry Parameters ==========
    int nx_, ny_;                      ///< Grid dimensions
    int totalCells_;                   ///< nx * ny
    double cellSizeX_, cellSizeY_;     ///< Cell side lengths
    int maxOccupancy_;                 ///< Maximum capacity per cell
    std::array<double, 2> boxSize_;    ///< Box boundaries (cached)
    
    // ========== Cell Data (SoA in cell) ==========
    /// Cell occupancy [totalCells]
    std::vector<int> occupancy_;
    
    /// Relative coordinates [totalCells][maxOccupancy] (flattened storage)
    std::vector<double> relPosX_;
    std::vector<double> relPosY_;
    
    /// Radius [totalCells][maxOccupancy] (key for polydisperse systems)
    std::vector<double> cellRadius_;
    
    /// Global particle ID [totalCells][maxOccupancy] (optional, for debugging)
    std::vector<int> particleID_;
    
    // ========== Neighbor Table (precomputed) ==========
    /// 3x3 neighbors [totalCells][9]
    /// Encoding: 0-8 corresponds to:
    ///   6  7  8
    ///   3  4  5
    ///   0  1  2
    std::vector<std::array<int, 9>> neighbors_;
    
    // ========== Helper Methods ==========
    /// Cell ID mapping: (ix, iy) -> cellID
    int cellIndex(int ix, int iy) const;
    
    /// Cell ID reverse mapping: cellID -> (ix, iy)
    std::pair<int, int> cellCoords(int cellID) const;
    
    /// Build neighbor table (considering PBC)
    void buildNeighborTable();
    
    /// Swap-tail deletion (O(1))
    void swapRemove(int cellID, int slotK);
    
    /// Automatically compute optimal grid size
    static std::pair<int, int> computeOptimalGrid(
        const std::array<double, 2>& boxSize, 
        double radiusMax, 
        int numParticles);
    
    /// 2D index flattening: cellID * maxOccupancy + slotK
    int flatIndex(int cellID, int slotK) const {
        return cellID * maxOccupancy_ + slotK;
    }
};

// ========== Template Method Implementation ==========

template<typename Func>
void CellGrid::forEachNeighbor(
    int activeCellID,
    int activeSlotK,
    const std::vector<int>& directions,
    Func callback) const
{
    // Get active particle's relative coordinates and radius
    int activeSlot = flatIndex(activeCellID, activeSlotK);
    double activeRelX = relPosX_[activeSlot];
    double activeRelY = relPosY_[activeSlot];
    
    // Get active particle's cell center coordinates
    std::pair<int, int> activeCoords = cellCoords(activeCellID);
    int activeIx = activeCoords.first;
    int activeIy = activeCoords.second;
    double activeCenterX = (activeIx + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
    double activeCenterY = (activeIy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
    
    // Active particle's global coordinates
    double activeGlobalX = activeCenterX + activeRelX;
    double activeGlobalY = activeCenterY + activeRelY;
    
    // Iterate over neighbors in specified directions
    for (int dir : directions) {
        int neighborCellID = neighbors_[activeCellID][dir];
        int ocp = occupancy_[neighborCellID];
        
        // Get neighbor cell's center coordinates
        std::pair<int, int> neighborCoords = cellCoords(neighborCellID);
        int neighborIx = neighborCoords.first;
        int neighborIy = neighborCoords.second;
        double neighborCenterX = (neighborIx + 0.5) * cellSizeX_ - boxSize_[0] / 2.0;
        double neighborCenterY = (neighborIy + 0.5) * cellSizeY_ - boxSize_[1] / 2.0;
        
        // Iterate over all particles in neighbor cell
        for (int k = 0; k < ocp; ++k) {
            int targetSlot = flatIndex(neighborCellID, k);
            double targetRelX = relPosX_[targetSlot];
            double targetRelY = relPosY_[targetSlot];
            double targetRadius = cellRadius_[targetSlot];
            
            // Target particle's global coordinates
            double targetGlobalX = neighborCenterX + targetRelX;
            double targetGlobalY = neighborCenterY + targetRelY;
            
            // Calculate relative displacement (PBC automatically handled)
            double dx = targetGlobalX - activeGlobalX;
            double dy = targetGlobalY - activeGlobalY;
            
            // PBC minimum image convention
            dx -= boxSize_[0] * std::round(dx / boxSize_[0]);
            dy -= boxSize_[1] * std::round(dy / boxSize_[1]);
            
            // Call callback function
            callback(neighborCellID, k, dx, dy, targetRadius);
        }
    }
}
