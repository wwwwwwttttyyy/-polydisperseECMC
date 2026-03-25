#pragma once
#include <cmath>
#include <limits>



inline double apply_pbc(double x, double L)
{
    return x - L * std::floor(x / L + 0.5);
}   
class IGeometry {
public:
    virtual ~IGeometry() = default;
    
    /**
     * @brief 计算事件距离（Event Chain Monte Carlo）
     * @param dx 目标粒子相对于活跃粒子的位置差（x2 - x1）
     * @param dy 目标粒子相对于活跃粒子的位置差（y2 - y1）
     * @param r1 活跃粒子半径
     * @param r2 目标粒子半径
     * @param vx, vy 事件链方向（通常是单位向量，如 (1,0) 表示 +x）
     * @return 事件距离 L（活跃粒子需移动的距离才碰撞）
     * 
     * 物理含义：
     *  - 活跃粒子从原点沿 (vx, vy) 移动距离 L
     *  - 目标粒子固定在 (dx, dy)
     *  - L 时刻两粒子刚好接触：|(L*vx, L*vy) - (dx, dy)| = r1 + r2
     * 
     * 返回值：
     *  - 正值：事件距离（未来会碰撞）
     *  - inf：不会碰撞（远离或平行）
     *  - 负值：不应该出现（说明已经重叠）
     */
    virtual double XeventDistance(
        double dx, double dy,
        double r1, double r2) const = 0;
    
    virtual double YeventDistance(
        double dx, double dy,
        double r1, double r2) const = 0;
};


/**
 * @brief 2D （Hard Disk）event distance
 * 
 * 碰撞条件：|r1(t) - r2| = r1 + r2
 * 求解二次方程：|r12 + v*t|² = (r1+r2)²
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
        
        // 快速剔除：Y 方向太远
        if (dy_sq >= sigma_sq) {
            return INF;
        }
        
        double lcoll = dx - std::sqrt(sigma_sq - dy_sq);
        return lcoll;
    }
    
    inline double YeventDistance(
        double dx, double dy,
        double r1, double r2) const override
    {
        double sigma = r1 + r2;
        double sigma_sq = sigma * sigma;
        double dx_sq = dx * dx;
        
        // 快速剔除：X 方向太远
        if (dx_sq >= sigma_sq) {
            return INF;
        }
        
        double lcoll = dy - std::sqrt(sigma_sq - dx_sq);
        return lcoll;
    }

    
private:
    static constexpr double INF = std::numeric_limits<double>::infinity();
    static constexpr double EPSILON = 1e-12;
};


/**
 * @brief 软球势能几何（未来扩展）
 * 
 * 势能：U(r) = ε * ((σ/r)^12 - (σ/r)^6)
 * 需要数值求解力为零的时刻
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
