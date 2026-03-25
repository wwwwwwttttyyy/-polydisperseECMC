#include "ecmc.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
#ifdef _WIN32
    #include <direct.h>
#else
    #include <sys/stat.h>
    #include <sys/types.h>
#endif

void makeDir(const char* name) {
#ifdef _WIN32
    _mkdir(name);   // Windows
#else
    mkdir(name, 0755);  // Linux / macOS
#endif
}

// 全局计数器，防止日志爆炸
EcmcEngine::EcmcEngine(System& sys, const AlgorithmPara& algo, const PhysicalPara& phys)
    : sys_(sys),
      grid_(sys),
      geom_(),
      rng_(algo.get_rng_seed()),
      algo_(algo),
      phys_(phys),
      pressure_observer_(),
      collision_count_(0)
{
    grid_.rebuild(sys_);
    displacement_max_ = calculateDisplacementMax();
    
    // 初始化观测器
    // pressure_observer_.
    pressure_observer_.initialize(
        sys_.size(),
        sys_.boxsize[0] * sys_.boxsize[1],
        sys_.radiusMean,
        algo_.get_chain_length(),
        algo_.get_sample_interval_pressure()
    );
    
    std::cout << "ECMC Engine: " << sys_.size() << " particles, "
              << algo_.get_n_steps() << " chains\n";
}

void EcmcEngine::run() {
    std::cout << "Running ECMC...\n";
    std::string outputDir = algo_.get_outputDir();
    makeDir(outputDir.c_str());
    std::string snapDir = outputDir + "/snapshot";
    makeDir(snapDir.c_str());

    std::string pressureDir = outputDir + "/pressure";
    makeDir(pressureDir.c_str());

    // 主循环：遍历所有事件链
     for (long long step = 0; step < algo_.get_n_steps(); ++step) {
        int direction = step % 2;
        std::pair<int, int> activeParticle = grid_.chooseActive(rng_);
        int cellID = activeParticle.first;
        int slotK = activeParticle.second;
        double chainLength = algo_.get_chain_length();
        double distanceToGo = chainLength;
        pressure_observer_.recordChain(direction);
        
        while (distanceToGo > 0) {

            double searchDistance = std::min(distanceToGo, displacement_max_);

            // 搜索碰撞
            CollisionEvent event;
            if (direction == 0) {
                event = findCollisionX(grid_, cellID, slotK, geom_, searchDistance);
            } else {
                event = findCollisionY(grid_, cellID, slotK, geom_, searchDistance);
            }

            // 决定实际移动距离
            double actualDisplacement = std::min({event.distance, distanceToGo, displacement_max_});
            actualDisplacement = std::max(0.0, actualDisplacement); // 非负
            int oldCellID = cellID;      
            int oldSlotK = slotK;
            std::pair<int, int> newPos;
            if (direction == 0) {
                newPos = grid_.moveParticle(sys_, cellID, slotK, actualDisplacement, 0.0);
            } else {
                newPos = grid_.moveParticle(sys_, cellID, slotK, 0.0, actualDisplacement);
            }
            int currentCellID = newPos.first;
            int currentSlotK = newPos.second;
            
            distanceToGo -= actualDisplacement;
            constexpr double EPS = 1e-12;
            bool is_collision = event.hasCollision() && (std::abs(actualDisplacement - event.distance) < EPS);

            if (is_collision) {

                if (currentCellID != oldCellID && event.targetCellID == oldCellID) {
                    int movedParticleOriginalSlot = grid_.getOccupancy(oldCellID);
                    if (event.targetSlotK == movedParticleOriginalSlot) {
                        event.targetSlotK = oldSlotK;
                    }
                }
                pressure_observer_.recordCollision(event.deltaIJ, direction);
                ++collision_count_;
                cellID = event.targetCellID;
                slotK = event.targetSlotK;
            } else {
                bool stopped_by_dmax = std::abs(actualDisplacement - displacement_max_) < EPS;
                cellID = currentCellID;
                slotK = currentSlotK;
            }
        }
        
        
        // 周期性压强采样
        if (algo_.get_sample_interval_pressure() > 0) {
            pressure_observer_.sampleAndReport();
        }
        
        // 压强弛豫完成后清除累积数据(在采样之后)
        if (algo_.get_n_equilibration_pressure() > 0 && step + 1 == algo_.get_n_equilibration_pressure()) {
            std::cout << "=== Pressure equilibration complete (" << (step + 1) << " chains) ===\n";
            std::cout << "Clearing pressure data, starting production sampling...\n";
            pressure_observer_.clear();
        }
        
        if (algo_.get_sample_interval_snapshot() > 0 &&
            step % algo_.get_sample_interval_snapshot() == 0)
        {
            char filename[256];
            sprintf(filename, "%s/snapshot_%09lld.dat", snapDir.c_str(), step);

            sys_.saveConfig(filename);
        }
        

    }
    // 末态重叠检查
    if (sys_.checkOverlap(1e-12)) {
        std::cerr << "Warning: Overlapping particles detected in final configuration!\n";
    }
    else {
        std::cout << "No overlaps detected in final configuration.\n";
    }
    std::cout << "Completed: " << collision_count_ << " collisions\n";
    



    // 输出压强统计数据
    pressure_observer_.writePressureData(pressureDir + "/pressure_data.txt");
}

double EcmcEngine::calculateDisplacementMax() const {
    double cellSizeX = grid_.getCellSizeX();
    double cellSizeY = grid_.getCellSizeY();
    double boxX = sys_.boxsize[0];
    double boxY = sys_.boxsize[1];
    double radiusMax = sys_.radiusMax;
    
    double limit = std::min({cellSizeX, cellSizeY, boxX / 2.0, boxY / 2.0});
    double displacement_max = limit - 2.0 * radiusMax;

    
    if (displacement_max <= 0) {
        throw std::runtime_error("displacement_max <= 0! System too dense or cell too small.");
    }
    
    return displacement_max;
}
