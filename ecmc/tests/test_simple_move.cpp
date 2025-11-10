#include "system.h"
#include "ecmc.h"
#include "para.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "========== Simple Movement Test ==========" << std::endl;
    std::cout << "Test: 1 particle, 1 chain, X direction, length 1.0" << std::endl;
    
    // 1. Create system
    System sys;
    sys.boxsize = {10.0, 10.0};
    sys.addParticle(2.0, 5.0, 0.5, 0);
    
    std::cout << "\nInitial position: x=" << sys.x[0] << ", y=" << sys.y[0] << std::endl;
    
    // 2. Configure parameters
    AlgorithmPara algo;
    algo.set_n_steps(1);            // Run only 1 chain
    algo.set_chain_length(1.0);     // Chain length 1.0
    algo.set_sample_interval_pressure(100);
    algo.set_sample_interval_snapshot(100);
    algo.set_rng_seed(12345);
    
    PhysicalPara phys;
    phys.set_n_atoms(1);
    phys.set_box({10.0, 10.0});
    
    // 3. Create ECMC engine
    EcmcEngine engine(sys, algo, phys);
    
    // 4. Run simulation
    std::cout << "\nRunning ECMC..." << std::endl;
    engine.run();
    
    // 5. Check results
    std::cout << "\nFinal position: x=" << sys.x[0] << ", y=" << sys.y[0] << std::endl;
    
    // Expected: step 0 direction = 0 (X direction)
    // Should move to x = 2.0 + 1.0 = 3.0, y = 5.0
    double expected_x = 3.0;
    double expected_y = 5.0;
    
    std::cout << "Expected position: x=" << expected_x << ", y=" << expected_y << std::endl;
    
    bool success = (std::abs(sys.x[0] - expected_x) < 1e-6 && 
                    std::abs(sys.y[0] - expected_y) < 1e-6);
    
    if (success) {
        std::cout << "\n✓ TEST PASSED" << std::endl;
        return 0;
    } else {
        std::cout << "\n❌ TEST FAILED" << std::endl;
        std::cout << "  dx = " << (sys.x[0] - expected_x) << std::endl;
        std::cout << "  dy = " << (sys.y[0] - expected_y) << std::endl;
        return 1;
    }
}
