#pragma once
#include <cmath>
#include <limits>

/**
 * @brief Geometry Collision Interface (Strategy Pattern)
 * 
 * Design goals:
 *  - Support different physical models (hard disks, soft spheres, ellipsoids, etc.)
 *  - Decouple collision mathematics from ECMC main logic
 *  - Easy testing and validation
 */
class IGeometry {
public:
    virtual ~IGeometry() = default;
    
    /**
     * @brief Calculate event distance (Event Chain Monte Carlo)
     * @param dx Target particle position difference relative to active (x2 - x1)
     * @param dy Target particle position difference relative to active (y2 - y1)
     * @param r1 Active particle radius
     * @param r2 Target particle radius
     * @param vx, vy Event chain direction (usually unit vector, e.g., (1,0) for +x)
     * @return Event distance L (distance active particle must move to collide)
     * 
     * Physical meaning:
     *  - Active particle moves from origin along (vx, vy) distance L
     *  - Target particle fixed at (dx, dy)
     *  - At time L, particles just touch: |(L*vx, L*vy) - (dx, dy)| = r1 + r2
     * 
     * Return values:
     *  - Positive: Event distance (will collide in future)
     *  - inf: Will not collide (moving away or parallel)
     *  - Negative: Should not occur (indicates overlap)
     */
    virtual double XeventDistance(
        double dx, double dy,
        double r1, double r2) const = 0;
    
    virtual double YeventDistance(
        double dx, double dy,
        double r1, double r2) const = 0;
};


/**
 * @brief 2D Hard Disk event distance
 * 
 * Collision condition: |r1(t) - r2| = r1 + r2
 * Solve quadratic equation: |r12 + v*t|^2 = (r1+r2)^2
 */
class HardDiskGeometry : public IGeometry {
public:
    inline double XeventDistance(
        double dx, double dy,
        double r1, double r2) const override
    {
        double sigma = r1 + r2;
        double sigma_sq = sigma * sigma;
        double dy_sq = dy * dy;
        
        // Fast rejection: Y direction too far
        if (dy_sq >= sigma_sq) {
            return INF;
        }
        
        double lcoll = dx - std::sqrt(sigma_sq - dy_sq);
        return lcoll > 0.0 ? lcoll : 0.0;
    }
    
    inline double YeventDistance(
        double dx, double dy,
        double r1, double r2) const override
    {
        double sigma = r1 + r2;
        double sigma_sq = sigma * sigma;
        double dx_sq = dx * dx;
        
        // Fast rejection: X direction too far
        if (dx_sq >= sigma_sq) {
            return INF;
        }
        
        double lcoll = dy - std::sqrt(sigma_sq - dx_sq);
        return lcoll > 0.0 ? lcoll : 0.0;
    }
    
private:
    static constexpr double INF = std::numeric_limits<double>::infinity();
    static constexpr double EPSILON = 1e-12;
};


/**
 * @brief Soft Sphere Potential Geometry (Future extension)
 * 
 * Potential: U(r) = epsilon * ((sigma/r)^12 - (sigma/r)^6)
 * Requires numerical solution for zero force moment
 */
class SoftSphereGeometry : public IGeometry {
public:
    SoftSphereGeometry(double epsilon, double sigma)
        : epsilon_(epsilon), sigma_(sigma) {}
    
    inline double XeventDistance(
        double dx, double dy,
        double r1, double r2) const override
    {
        return std::numeric_limits<double>::infinity();
    }
    
    inline double YeventDistance(
        double dx, double dy,
        double r1, double r2) const override
    {
        return std::numeric_limits<double>::infinity();
    }
    
private:
    double epsilon_;
    double sigma_;
};
