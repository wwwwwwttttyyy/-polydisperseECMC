#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>

class System {
public:
    // ========== Construction and Initialization ==========
    System() = default;
    explicit System(size_t n, const std::array<double, 2>& box_);  ///< Preallocate n particles and set box boundaries
    
    std::vector<double> x;        ///< x coordinates [N]
    std::vector<double> y;        ///< y coordinates [N]
    std::vector<double> radius;   ///< radius [N] (key for polydisperse systems)
    std::vector<int> typeID;      ///< type ID [N] (optional, for statistics)
    std::array<double, 2> boxsize;    ///< {Lx, Ly}, should be set externally based on PhysicalPara, assumes box center at origin in simulation
    // Store max radius for cell list construction
    double radiusMax = 0.0;       ///< max(radius)
    double radiusMean = 0.0;      ///< mean(radius)
    double radiusMin = std::numeric_limits<double>::infinity();  ///< min(radius)

    size_t size() const { return x.size(); }
    
    // Preallocate memory
    void reserve(size_t n);

    // Clear all data and free space
    void clear();
    
    // ========== Particle Operations ==========
    void addParticle(double x_, double y_, double r_, int tid = 0);
    void removeParticle(size_t index);
    void updateMaxMeanRadius();
    
    // ========== File I/O ==========
    /**
     * @brief Save system configuration to text file (human readable)
     * @param filename file name
     * 
     * Format:
     *   Line 1: Lx Ly
     *   Following lines: x y radius
     */
    void saveConfig(const std::string& filename) const;
    
    /**
     * @brief Load system configuration from text file
     * @param filename file name
     * 
     * Same format as saveConfig
     * Clears current system and loads new data
     */
    void loadConfig(const std::string& filename);
    
    /**
     * @brief Save system configuration to binary file (efficient)
     * @param filename file name
     * 
     * Format (all data as double):
     *   - Lx, Ly (box dimensions)
     *   - N (number of particles, stored as double)
     *   - x[0], y[0], r[0]
     *   - x[1], y[1], r[1]
     *   - ...
     */
    void saveConfigBinary(const std::string& filename) const;
    
    /**
     * @brief Load system configuration from binary file
     * @param filename file name
     */
    void loadConfigBinary(const std::string& filename);
    
    // Box operations and global PBC implementation
    // Initial configuration overlap check. false: no overlap, true: overlap detected
    bool checkOverlap() const;
    // Wrap back to main box
    void applyPBC(size_t i);
    // Global PBC
    std::array<double, 2> displacement(size_t i, size_t j) const;
    // Generalized distance
    double distance(size_t i, size_t j) const;
    
    // ========== Statistical Calculations ==========
    /**
     * @brief Calculate system packing fraction
     * @return φ = Σ(πr²) / (Lx * Ly)
     */
    double getPackingFraction() const;
};
