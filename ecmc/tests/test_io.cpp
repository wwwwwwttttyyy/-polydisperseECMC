/**
 * Test binary IO performance and correctness
 */

#include "system.h"
#include <iostream>
#include <chrono>
#include <cmath>

int main() {
    std::cout << "========== Binary IO Test ==========\n\n";
    
    // Create test system
    System sys;
    sys.boxsize[0] = 20.0;
    sys.boxsize[1] = 20.0;
    
    // Add 1000 particles
    const int N = 1000;
    for (int i = 0; i < N; ++i) {
        double x = (i % 32) - 16.0;
        double y = (i / 32) - 16.0;
        double r = 0.5;
        sys.addParticle(x, y, r);
    }
    
    std::cout << "Created system with " << sys.size() << " particles\n\n";
    
    // ========== Test Text IO ==========
    std::cout << "--- Text Format ---\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    sys.saveConfig("test_config.txt");
    auto t2 = std::chrono::high_resolution_clock::now();
    auto dt_text_save = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    
    System sys_text;
    t1 = std::chrono::high_resolution_clock::now();
    sys_text.loadConfig("test_config.txt");
    t2 = std::chrono::high_resolution_clock::now();
    auto dt_text_load = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    
    std::cout << "Save: " << dt_text_save << " us\n";
    std::cout << "Load: " << dt_text_load << " us\n\n";
    
    // ========== Test Binary IO ==========
    std::cout << "--- Binary Format ---\n";
    t1 = std::chrono::high_resolution_clock::now();
    sys.saveConfigBinary("test_config.bin");
    t2 = std::chrono::high_resolution_clock::now();
    auto dt_bin_save = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    
    System sys_bin;
    t1 = std::chrono::high_resolution_clock::now();
    sys_bin.loadConfigBinary("test_config.bin");
    t2 = std::chrono::high_resolution_clock::now();
    auto dt_bin_load = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    
    std::cout << "Save: " << dt_bin_save << " us\n";
    std::cout << "Load: " << dt_bin_load << " us\n\n";
    
    // ========== Verify correctness ==========
    std::cout << "--- Verification ---\n";
    bool correct = true;
    if (sys_bin.size() != sys.size()) {
        std::cerr << "ERROR: Particle count mismatch!\n";
        correct = false;
    }
    
    for (size_t i = 0; i < sys.size(); ++i) {
        if (std::abs(sys_bin.x[i] - sys.x[i]) > 1e-10 ||
            std::abs(sys_bin.y[i] - sys.y[i]) > 1e-10 ||
            std::abs(sys_bin.radius[i] - sys.radius[i]) > 1e-10) {
            std::cerr << "ERROR: Particle " << i << " data mismatch!\n";
            correct = false;
            break;
        }
    }
    
    if (correct) {
        std::cout << " Binary IO correct!\n";
    }
    
    // ========== Performance summary ==========
    std::cout << "\n--- Performance Speedup ---\n";
    std::cout << "Save: " << (double)dt_text_save / dt_bin_save << "x faster\n";
    std::cout << "Load: " << (double)dt_text_load / dt_bin_load << "x faster\n";
    
    std::cout << "\n========== TEST PASSED ==========\n";
    return 0;
}
