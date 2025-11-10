#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>

class System {
public:
    // ========== Construction & Initialization ==========
    System() = default;
    explicit System(size_t n, const std::array<double, 2>& box_);  ///< Pre-allocate n particles
    
    std::vector<double> x;        ///< X coordinates [N]
    std::vector<double> y;        ///< Y coordinates [N]
    std::vector<double> radius;   ///< Radii [N] (polydispersity support)
    std::vector<int> typeID;      ///< Type ID [N] (optional)
    std::array<double, 2> boxsize;    ///< {Lx, Ly}, box centered at origin
    // Store max/mean/min radius for cell grid construction
    double radiusMax = 0.0;       ///< max(radius)
    double radiusMean = 0.0;      ///< mean(radius)
    double radiusMin = std::numeric_limits<double>::infinity();  ///< min(radius)

    size_t size() const { return x.size(); }
    
    void reserve(size_t n);
    void clear();
    
    // ========== Particle Operations ==========
    void addParticle(double x_, double y_, double r_, int tid = 0);
    void removeParticle(size_t index);
    void updateMaxMeanRadius();
    
    // ========== File I/O ==========
    /**
     * @brief Save configuration to text file (human-readable)
     * Format: Line 1: Lx Ly; Following lines: x y radius
     */
    void saveConfig(const std::string& filename) const;
    
    /**
     * @brief Load configuration from text file
     */
    void loadConfig(const std::string& filename);
    
    /**
     * @brief Save configuration to binary file (efficient)
     * Format: Lx, Ly, N, then x[i] y[i] r[i] for each particle
     */
    void saveConfigBinary(const std::string& filename) const;
    
    /**
     * @brief Load configuration from binary file
     */
    void loadConfigBinary(const std::string& filename);
    
    // ========== Box Operations & PBC ==========
    bool checkOverlap() const;  ///< Check for overlaps, return true if any found
    void applyPBC(size_t i);    ///< Apply PBC to particle i
    std::array<double, 2> displacement(size_t i, size_t j) const;  ///< Displacement with PBC
    double distance(size_t i, size_t j) const;  ///< Distance with PBC
};