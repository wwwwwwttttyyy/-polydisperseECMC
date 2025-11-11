#pragma once
#include "system.h"
#include "cell.h"
#include "geometry.h"
#include "collision.h"
#include "para.h"
#include "pressureObserver.h"
#include <random>

/**
 * @brief ECMC (Event Chain Monte Carlo) main engine
 * 
 * Only responsible for main loop orchestration, all specific calculations done by components
 */
class EcmcEngine {
public:
    EcmcEngine(System& sys, const AlgorithmPara& algo, const PhysicalPara& phys);
    
    void run();
    
    long long getCollisionCount() const { return collision_count_; }
    PressureObserver& getPressureObserver() { return pressure_observer_; }
    const PressureObserver& getPressureObserver() const { return pressure_observer_; }
    const CellGrid& getGrid() const { return grid_; }
    
private:
    double calculateDisplacementMax() const;
    
    System& sys_;
    CellGrid grid_;
    HardDiskGeometry geom_;
    std::mt19937 rng_;
    const AlgorithmPara& algo_;
    const PhysicalPara& phys_;
    PressureObserver pressure_observer_;
    long long collision_count_;
    double displacement_max_;
};
