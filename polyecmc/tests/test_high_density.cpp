#include "system.h"
#include "ecmc.h"
#include "para.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>

// Check if two particles overlap
bool checkOverlap(const System& sys, size_t i, size_t j) {
    double dx = sys.x[j] - sys.x[i];
    double dy = sys.y[j] - sys.y[i];
    dx -= sys.boxsize[0] * std::round(dx / sys.boxsize[0]);
    dy -= sys.boxsize[1] * std::round(dy / sys.boxsize[1]);
    double dist = std::sqrt(dx*dx + dy*dy);
    double sigma = sys.radius[i] + sys.radius[j];
    return dist < (sigma - 1e-10);
}

// Check if the entire system has overlaps
bool checkAllNoOverlap(const System& sys) {
    size_t N = sys.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            if (checkOverlap(sys, i, j)) {
                std::cout << "Overlap detected between particles " << i 
                          << " and " << j << std::endl;
                double dx = sys.x[j] - sys.x[i];
                double dy = sys.y[j] - sys.y[i];
                dx -= sys.boxsize[0] * std::round(dx / sys.boxsize[0]);
                dy -= sys.boxsize[1] * std::round(dy / sys.boxsize[1]);
                double dist = std::sqrt(dx*dx + dy*dy);
                double sigma = sys.radius[i] + sys.radius[j];
                std::cout << "  Distance: " << dist << ", Sigma: " << sigma << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Output configuration to file (for Python visualization)
void saveConfiguration(const System& sys, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    // Write box size
    file << "# Box: " << sys.boxsize[0] << " " << sys.boxsize[1] << "\n";
    
    // Write particle count
    file << "# N: " << sys.size() << "\n";
    
    // Write particle positions and radii (x, y, radius)
    file << "# x y radius\n";
    for (size_t i = 0; i < sys.size(); ++i) {
        file << std::setprecision(10) << sys.x[i] << " " 
             << sys.y[i] << " " << sys.radius[i] << "\n";
    }
    
    file.close();
    std::cout << "Configuration saved to " << filename << std::endl;
}

// Generate square lattice initial configuration
void generateSquareLattice(System& sys, double packing_fraction, double radius) {
    sys.clear();
    
    // For square lattice, each particle occupies area a²
    // packing_fraction = π * r² / a²
    // Therefore a = r * sqrt(π / packing_fraction)
    
    double a = radius * std::sqrt(M_PI / packing_fraction);  // Lattice constant
    
    // First determine a reasonable box size
    int nx = 16;  // Number of particles in X direction
    int ny = 16;  // Number of particles in Y direction
    
    double Lx = nx * a;
    double Ly = ny * a;
    
    sys.boxsize = {Lx, Ly};
    
    std::cout << "Lattice parameters:" << std::endl;
    std::cout << "  Lattice constant a = " << a << " (should be > 2*r = " << 2*radius << ")" << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << std::endl;
    std::cout << "  Box size: [" << Lx << ", " << Ly << "]" << std::endl;
    
    // Generate square lattice (centered)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = (i + 0.5) * a - Lx / 2.0;
            double y = (j + 0.5) * a - Ly / 2.0;
            sys.addParticle(x, y, radius, 0);
        }
    }
    
    // Verify packing fraction
    double area = Lx * Ly;
    double particle_area = sys.size() * M_PI * radius * radius;
    double actual_packing = particle_area / area;
    
    std::cout << "  N = " << sys.size() << " particles" << std::endl;
    std::cout << "  Target packing fraction = " << packing_fraction << std::endl;
    std::cout << "  Actual packing fraction = " << actual_packing << std::endl;
}

// Generate random non-overlapping initial configuration (simple Monte Carlo method)
void generateRandomConfiguration(System& sys, int N, double radius, 
                                  double Lx, double Ly, int max_attempts = 10000) {
    std::mt19937 rng(42);  // Fixed random seed for reproducibility
    std::uniform_real_distribution<double> distX(-Lx/2.0, Lx/2.0);
    std::uniform_real_distribution<double> distY(-Ly/2.0, Ly/2.0);
    
    sys.clear();
    sys.boxsize = {Lx, Ly};
    
    int attempts = 0;
    while (sys.size() < static_cast<size_t>(N) && attempts < max_attempts) {
        double x = distX(rng);
        double y = distY(rng);
        
        // Check if overlaps with existing particles
        bool overlap = false;
        for (size_t i = 0; i < sys.size(); ++i) {
            double dx = x - sys.x[i];
            double dy = y - sys.y[i];
            dx -= Lx * std::round(dx / Lx);
            dy -= Ly * std::round(dy / Ly);
            double dist = std::sqrt(dx*dx + dy*dy);
            if (dist < 2.0 * radius + 0.01) {  // Leave some gap
                overlap = true;
                break;
            }
        }
        
        if (!overlap) {
            sys.addParticle(x, y, radius, 0);
            attempts = 0;  // Reset attempt count
        } else {
            attempts++;
        }
    }
    
    if (sys.size() < static_cast<size_t>(N)) {
        std::cout << "Warning: Only placed " << sys.size() << " out of " 
                  << N << " particles (system too dense)" << std::endl;
    }
}

int main() {
    std::cout << "========== High Density ECMC Test ==========" << std::endl;
    
    // Parameter settings
    double packing_fraction = 0.7;  // Target packing fraction
    double radius = 0.5;             // Particle radius
    
    std::cout << "System parameters:" << std::endl;
    std::cout << "  Target packing fraction = " << packing_fraction << std::endl;
    std::cout << "  Radius = " << radius << std::endl;
    
    // Generate square lattice initial configuration
    std::cout << "\nGenerating square lattice initial configuration..." << std::endl;
    System sys;
    generateSquareLattice(sys, packing_fraction, radius);
    
    // Modify to 32x32 grid (1024 particles)
    sys.clear();
    
    double a = radius * std::sqrt(M_PI / packing_fraction);
    int nx = 32;  // 32 x 32 = 1024 particles
    int ny = 32;
    
    double Lx = nx * a;
    double Ly = ny * a;
    sys.boxsize = {Lx, Ly};
    
    std::cout << "Lattice parameters:" << std::endl;
    std::cout << "  Lattice constant a = " << a << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << std::endl;
    std::cout << "  Box size: [" << Lx << ", " << Ly << "]" << std::endl;
    
    // Generate square lattice
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = (i + 0.5) * a - Lx / 2.0;
            double y = (j + 0.5) * a - Ly / 2.0;
            sys.addParticle(x, y, radius, 0);
        }
    }
    
    std::cout << "  N = " << sys.size() << " particles" << std::endl;
    double area = Lx * Ly;
    double particle_area = sys.size() * M_PI * radius * radius;
    double actual_packing = particle_area / area;
    std::cout << "  Target packing fraction = " << packing_fraction << std::endl;
    std::cout << "  Actual packing fraction = " << actual_packing << std::endl;
    
    int N = sys.size();
    
    // Check initial configuration
    if (!checkAllNoOverlap(sys)) {
        std::cout << "❌ Initial configuration has overlaps!" << std::endl;
        return 1;
    }
    std::cout << "\n✓ Initial configuration: No overlaps" << std::endl;
    
    // Save initial configuration
    saveConfiguration(sys, "config_initial.dat");
    
    // Configure ECMC parameters
    AlgorithmPara algo;
    algo.set_n_steps(1000000);       // Run 1 million chains (performance test)
    algo.set_chain_length(1.0);      // Chain length 1.0
    algo.set_sample_interval_pressure(100);
    algo.set_sample_interval_snapshot(10000);
    algo.set_rng_seed(12345);
    
    PhysicalPara phys;
    phys.set_n_atoms(N);
    phys.set_box({sys.boxsize[0], sys.boxsize[1]});
    
    // Run ECMC
    std::cout << "\n========== Running ECMC ==========" << std::endl;
    EcmcEngine engine(sys, algo, phys);
    
    engine.run();
    
    // Check final configuration
    std::cout << "\n========== Validation ==========" << std::endl;
    if (!checkAllNoOverlap(sys)) {
        std::cout << "❌ FAILED: Final configuration has overlaps!" << std::endl;
        saveConfiguration(sys, "config_final_FAILED.dat");
        return 1;
    }
    
    std::cout << "✓ Final configuration: No overlaps" << std::endl;
    
    // Save final configuration
    saveConfiguration(sys, "config_final.dat");
    
    std::cout << "\n========== TEST PASSED ==========" << std::endl;
    std::cout << "Initial and final configurations saved." << std::endl;
    std::cout << "Use Python script to visualize." << std::endl;
    
    return 0;
}
