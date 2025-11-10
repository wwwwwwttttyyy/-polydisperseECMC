#include "../include/system.h"
#include "../include/cell.h"
#include "../include/ecmc.h"
#include "../include/para.h"
#include <iostream>
#include <iomanip>
#include <cmath>

// Check if two particles overlap
bool checkOverlap(const System& sys, int i, int j) {
    double dx = sys.x[i] - sys.x[j];
    double dy = sys.y[i] - sys.y[j];
    
    // PBC
    if (dx > sys.boxsize[0] / 2) dx -= sys.boxsize[0];
    if (dx < -sys.boxsize[0] / 2) dx += sys.boxsize[0];
    if (dy > sys.boxsize[1] / 2) dy -= sys.boxsize[1];
    if (dy < -sys.boxsize[1] / 2) dy += sys.boxsize[1];
    
    double dist = std::sqrt(dx * dx + dy * dy);
    double sigma = sys.radius[i] + sys.radius[j];
    
    return dist < sigma - 1e-10;  // Tolerate floating point errors
}

// Check if all particles have no overlaps
bool checkAllNoOverlap(const System& sys) {
    for (size_t i = 0; i < sys.size(); ++i) {
        for (size_t j = i + 1; j < sys.size(); ++j) {
            if (checkOverlap(sys, i, j)) {
                std::cout << "ERROR: Particles " << i << " and " << j << " overlap!" << std::endl;
                
                double dx = sys.x[i] - sys.x[j];
                double dy = sys.y[i] - sys.y[j];
                
                // PBC
                if (dx > sys.boxsize[0] / 2) dx -= sys.boxsize[0];
                if (dx < -sys.boxsize[0] / 2) dx += sys.boxsize[0];
                if (dy > sys.boxsize[1] / 2) dy -= sys.boxsize[1];
                if (dy < -sys.boxsize[1] / 2) dy += sys.boxsize[1];
                
                double dist = std::sqrt(dx * dx + dy * dy);
                double sigma = sys.radius[i] + sys.radius[j];
                
                std::cout << "  Distance: " << dist << ", Sigma: " << sigma 
                          << ", Overlap: " << (sigma - dist) << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Output particle positions
void printPositions(const System& sys, const std::string& label) {
    std::cout << "\n" << label << ":\n";
    for (size_t i = 0; i < sys.size(); ++i) {
        std::cout << "  Particle " << i << ": ("
                  << std::fixed << std::setprecision(6)
                  << sys.x[i] << ", "
                  << sys.y[i] << "), r="
                  << sys.radius[i] << std::endl;
    }
}

int main() {
    std::cout << "========== ECMC Movement Test ==========" << std::endl;
    
    // 1. Create simple system: 2 particles
    System sys;
    sys.boxsize = {10.0, 10.0};
    
    // Particle 0: left
    sys.addParticle(2.0, 5.0, 0.5, 0);
    
    // Particle 1: right, distance just > σ
    sys.addParticle(4.2, 5.0, 0.5, 1);
    
    std::cout << "Box size: [" << sys.boxsize[0] << ", " << sys.boxsize[1] << "]" << std::endl;
    std::cout << "Number of particles: " << sys.size() << std::endl;
    
    printPositions(sys, "Initial positions");
    
    // Check initial state
    if (!checkAllNoOverlap(sys)) {
        std::cout << "\nERROR: Initial configuration has overlaps!" << std::endl;
        return 1;
    }
    std::cout << "\n✓ Initial configuration: No overlaps" << std::endl;
    
    // 2. Set algorithm parameters
    AlgorithmPara algo;
    algo.set_n_steps(10);                    // Run only 10 chains
    algo.set_chain_length(1.0);              // Short chain length
    algo.set_sample_interval_pressure(100);  // No sampling for now
    algo.set_sample_interval_snapshot(100);
    algo.set_rng_seed(12345);
    
    PhysicalPara phys;
    phys.set_n_atoms(2);
    phys.set_box({10.0, 10.0});
    
    // 3. Create ECMC engine
    EcmcEngine engine(sys, algo, phys);
    
    // Debug: Check System state
    std::cout << "\n========== DEBUG: System State ==========" << std::endl;
    std::cout << "sys.size() = " << sys.size() << std::endl;
    std::cout << "sys.radiusMax = " << sys.radiusMax << std::endl;
    std::cout << "sys.radiusMean = " << sys.radiusMean << std::endl;
    std::cout << "Particle 0: x=" << sys.x[0] << ", y=" << sys.y[0] << ", r=" << sys.radius[0] << std::endl;
    std::cout << "Particle 1: x=" << sys.x[1] << ", y=" << sys.y[1] << ", r=" << sys.radius[1] << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // 4. Run simulation
    std::cout << "\n========== Running ECMC ==========" << std::endl;
    engine.run();
    
    // 5. Check final state
    printPositions(sys, "Final positions");
    
    std::cout << "\n========== Validation ==========" << std::endl;
    if (!checkAllNoOverlap(sys)) {
        std::cout << "\nFAILED: Particles overlap after simulation!" << std::endl;
        return 1;
    }
    
    std::cout << "\nAll particles have valid positions (no overlaps)" << std::endl;
    
    // 6. Check if particles moved
    bool moved = false;
    double initial_pos = 2.0;  // Particle 0 initial position
    if (std::abs(sys.x[0] - initial_pos) > 1e-6) {
        moved = true;
    }
    
    if (moved) {
        std::cout << " Particles moved during simulation" << std::endl;
    } else {
        std::cout << "Warning: Particles did not move (may be expected for short runs)" << std::endl;
    }
    
    // 7. Check collision count
    long long collisions = engine.getCollisionCount();
    std::cout << "✓ Total collisions: " << collisions << std::endl;
    
    std::cout << "\n========== TEST PASSED ==========" << std::endl;
    return 0;
}
 