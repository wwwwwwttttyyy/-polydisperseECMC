#include "ecmc.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>

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
    
    // Initialize observers
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
    
    // Main loop: iterate over all event chains
    for (long long step = 0; step < algo_.get_n_steps(); ++step) {
        // Alternate direction (0=X, 1=Y)
        int direction = step % 2;
        
        // Choose active particle
        std::pair<int, int> activeParticle = grid_.chooseActive(rng_);
        int cellID = activeParticle.first;
        int slotK = activeParticle.second;
        
        // ========== Get event chain length ==========
        double chainLength = algo_.get_chain_length();
        double distanceToGo = chainLength;
        
        // Record chain count to PressureObserver
        pressure_observer_.recordChain(direction);
        
        // Collision loop: move along specified direction until chain ends
        while (distanceToGo > 0) {
            // Calculate current step search distance (limited by displacement_max)
            double searchDistance = std::min(distanceToGo, displacement_max_);
            
            // Search for collision
            CollisionEvent event;
            if (direction == 0) {
                event = findCollisionX(grid_, cellID, slotK, geom_, searchDistance);
            } else {
                event = findCollisionY(grid_, cellID, slotK, geom_, searchDistance);
            }
            
            // Ensure event distance is non-negative (numerical error protection)
            double eventDist = std::max(event.distance, 0.0);
            
            // Calculate actual displacement
            double actualDisplacement = std::min(std::min(eventDist, distanceToGo), displacement_max_);
            
            // Move particle: CellGrid will simultaneously update System global coordinates (with PBC)
            if (direction == 0) {
                // X direction movement
                std::pair<int, int> newPos = grid_.moveParticle(sys_, cellID, slotK, actualDisplacement, 0.0);
                cellID = newPos.first;
                slotK = newPos.second;
            } else {
                // Y direction movement
                std::pair<int, int> newPos = grid_.moveParticle(sys_, cellID, slotK, 0.0, actualDisplacement);
                cellID = newPos.first;
                slotK = newPos.second;
            }
            
            // Determine three cases, here are three physical quantities: lcoll, lremain, lmax
            // lcoll smallest: normal chain propagation;
            // lremain smallest: move lremain then chain ends
            // lmax smallest: approaching cell limit, move lmax then same particle continues

            // lmax smallest
            if (displacement_max_ < std::min(eventDist, distanceToGo)) 
            {
                distanceToGo -= displacement_max_;
                continue;

            }
            // lcoll smallest, normal chain propagation
            else if (eventDist < distanceToGo) {

                // Record collision to PressureObserver (deltaIJ already calculated in findCollision)
                pressure_observer_.recordCollision(event.deltaIJ, direction);
                ++collision_count_;
                
                // Switch active particle
                cellID = event.targetCellID;
                slotK = event.targetSlotK;
                
                distanceToGo -= eventDist;
                continue;
            }
            else 
            // lremain smallest
            {
                break;
            }
        }
        
        // Progress output
        if ((step + 1) % 10000 == 0) {
            std::cout << "  " << (step + 1) << " / " << algo_.get_n_steps()
                      << " chains (" << collision_count_ << " collisions)\n";
        }
        
        // Periodic pressure sampling
        if ((step + 1) % algo_.get_sample_interval_pressure() == 0) {
            pressure_observer_.sampleAndReport();
        }
        
        // Clear accumulated data after pressure equilibration (after sampling)
        if (algo_.get_n_equilibration_pressure() > 0 && step + 1 == algo_.get_n_equilibration_pressure()) {
            std::cout << "=== Pressure equilibration complete (" << (step + 1) << " chains) ===\n";
            std::cout << "Clearing pressure data, starting production sampling...\n";
            pressure_observer_.clear();
        }
        
        // Periodic configuration snapshot output
        if (algo_.get_sample_interval_snapshot() > 0 && 
            (step + 1) % algo_.get_sample_interval_snapshot() == 0) {
            // Zero-padded filename: snapshot_00000100000.dat
            char filename[64];
            std::sprintf(filename, "snapshot_%011lld.dat", step + 1);
            sys_.saveConfigBinary(filename);
            std::cout << "Snapshot saved: " << filename << "\n";
        }
        
        // Periodic overlap check
        if ((step + 1) % 100000 == 0) {
            if (sys_.checkOverlap()) {
                throw std::runtime_error("Overlap detected at step " + std::to_string(step + 1));
            }
        }
    }
    
    std::cout << "Completed: " << collision_count_ << " collisions\n";
    
    // Output pressure statistics data
    pressure_observer_.writePressureData("pressure_data.txt");
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
