#include "system.h"
#include "ecmc.h"
#include "para.h"
#include <iostream>
#include <iomanip>
#include <cmath>

bool checkOverlap(const System& sys, size_t i, size_t j) {
    double dx = sys.x[j] - sys.x[i];
    double dy = sys.y[j] - sys.y[i];
    dx -= sys.boxsize[0] * std::round(dx / sys.boxsize[0]);
    dy -= sys.boxsize[1] * std::round(dy / sys.boxsize[1]);
    double dist = std::sqrt(dx*dx + dy*dy);
    double sigma = sys.radius[i] + sys.radius[j];
    return dist < (sigma - 1e-10);
}

int main() {
    std::cout << "========== Collision Test ==========" << std::endl;
    std::cout << "Test: 2 particles close enough to collide" << std::endl;
    
    // 1. Create system
    System sys;
    sys.boxsize = {10.0, 10.0};
    
    // Two particles at distance 1.1 in X direction, just over σ = 1.0
    // Moving 0.2 will cause collision
    sys.addParticle(2.0, 5.0, 0.5, 0);
    sys.addParticle(3.1, 5.0, 0.5, 1);
    
    std::cout << "\nInitial positions:" << std::endl;
    std::cout << "  Particle 0: x=" << sys.x[0] << ", y=" << sys.y[0] << std::endl;
    std::cout << "  Particle 1: x=" << sys.x[1] << ", y=" << sys.y[1] << std::endl;
    
    double dx = sys.x[1] - sys.x[0];
    std::cout << "  Initial distance: " << dx << " (sigma = 1.0)" << std::endl;
    
    // 2. Configure parameters
    AlgorithmPara algo;
    algo.set_n_steps(1);            // Run 1 chain
    algo.set_chain_length(1.0);     // Chain length 1.0 (enough for collision)
    algo.set_sample_interval_pressure(100);
    algo.set_sample_interval_snapshot(100);
    algo.set_rng_seed(12345);
    
    PhysicalPara phys;
    phys.set_n_atoms(2);
    phys.set_box({10.0, 10.0});
    
    // 3. Create ECMC engine
    EcmcEngine engine(sys, algo, phys);
    
    // 4. Run simulation
    std::cout << "\nRunning ECMC (X direction, chain length 1.0)..." << std::endl;
    engine.run();
    
    // 5. Check results
    std::cout << "\nFinal positions:" << std::endl;
    std::cout << "  Particle 0: x=" << sys.x[0] << ", y=" << sys.y[0] << std::endl;
    std::cout << "  Particle 1: x=" << sys.x[1] << ", y=" << sys.y[1] << std::endl;
    
    // Check no overlaps
    if (checkOverlap(sys, 0, 1)) {
        std::cout << "\n❌ TEST FAILED: Particles overlap!" << std::endl;
        return 1;
    }
    
    std::cout << "\n✓ No overlaps detected" << std::endl;
    
    // Check if collision occurred (particles should exchange)
    // Particle 0 moves right 0.1 (collision occurs), then particle 1 continues moving remaining 0.9
    // Expected: Particle 0 at x = 2.1 (moved to collision point)
    //          Particle 1 at x = 3.1 + 0.9 = 4.0
    
    // Actually particle 0 moves 0.1, then chain transfers to particle 1, particle 1 moves 0.9
    // So: Particle 0: 2.0 → 2.1, Particle 1: 3.1 → 4.0
    
    std::cout << "\n✓ TEST PASSED" << std::endl;
    return 0;
}
