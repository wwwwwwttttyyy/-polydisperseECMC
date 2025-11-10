/**
 * ECMC (Event-Chain Monte Carlo) Demonstration
 * Example: 256 hard disks, square lattice to liquid equilibration
 */

#include "system.h"
#include "cell.h"
#include "geometry.h"
#include "collision.h"
#include "ecmc.h"
#include "para.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

void generateSquareLattice(System& sys, int nx, int ny, double phi, double polydispersity = 0.0, unsigned int seed = 123) {
    int N = nx * ny;
    double r_mean = 1.0;
    double a = std::sqrt(M_PI * r_mean * r_mean / phi);
    double Lx = nx * a, Ly = ny * a;
    
    sys.x.resize(N);
    sys.y.resize(N);
    sys.radius.resize(N);
    sys.boxsize[0] = Lx;
    sys.boxsize[1] = Ly;
    
    // Generate polydisperse radii
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(r_mean * (1.0 - polydispersity), 
                                                 r_mean * (1.0 + polydispersity));
    
    // Perfect square lattice (no perturbation for high-density systems)
    int idx = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            sys.x[idx] = (i + 0.5) * a - Lx / 2.0;
            sys.y[idx] = (j + 0.5) * a - Ly / 2.0;
            sys.radius[idx] = dist(rng);
            ++idx;
        }
    }
    
    sys.updateMaxMeanRadius();
    
    std::cout << "Square lattice: " << nx << "x" << ny << " = " << N << " particles\n";
    std::cout << "  Box: " << Lx << " x " << Ly << ", phi = " << phi << "\n";
    if (polydispersity > 0) {
        std::cout << "  Polydispersity: +/-" << (polydispersity * 100) << "% "
                  << "(r in [" << r_mean * (1.0 - polydispersity) << ", " 
                  << r_mean * (1.0 + polydispersity) << "])\n";
        std::cout << "  r_mean = " << sys.radiusMean << ", r_max = " << sys.radiusMax 
                  << ", r_min = " << sys.radiusMin << "\n";
    }
    std::cout << "\n";
}

void saveConfig(const System& sys, const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: Cannot open " << filename << "\n";
        return;
    }
    ofs << sys.boxsize[0] << " " << sys.boxsize[1] << "\n";
    for (size_t i = 0; i < sys.x.size(); ++i) {
        ofs << sys.x[i] << " " << sys.y[i] << " " << sys.radius[i] << "\n";
    }
}

int main() {
    std::cout << "========== ECMC Demonstration ==========\n\n";
    
    // Configuration
    const int nx = 16, ny = 16;
    const double phi = 0.7;
    const int n_chains = 500000;
    const unsigned int rng_seed = 42;
    const double polydispersity = 0.02;  // Weak polydispersity: +/- 2% (radius 0.98-1.02)
    
    // Initialize system
    System sys;
    generateSquareLattice(sys, nx, ny, phi, polydispersity, rng_seed);
    
    if (sys.checkOverlap()) {
        std::cerr << "Error: Initial overlaps detected!\n";
        return 1;
    }
    std::cout << "Initial: No overlaps\n";
    saveConfig(sys, "config_initial.dat");
    std::cout << "Saved to config_initial.dat\n\n";
    
    // Setup ECMC
    std::cout << "========== Running ECMC ==========\n";
    
    AlgorithmPara algo_para;
    algo_para.set_n_steps(n_chains);
    algo_para.set_chain_length(1.0);
    algo_para.set_sample_interval_pressure(100000);
    algo_para.set_sample_interval_snapshot(n_chains);
    algo_para.set_rng_seed(rng_seed);
    algo_para.set_n_equilibration_pressure(100000);
    
    PhysicalPara phys_para;
    phys_para.set_n_atoms(nx * ny);
    phys_para.set_VolumeFraction(phi);
    phys_para.set_box({sys.boxsize[0], sys.boxsize[1]});
    
    EcmcEngine algo(sys, algo_para, phys_para);
    
    // Print Cell Grid info
    std::cout << "Cell Grid: " << algo.getGrid().getNx() << " x " 
              << algo.getGrid().getNy() << " = " 
              << algo.getGrid().getTotalCells() << " cells\n";
    std::cout << "  Cell size: " << algo.getGrid().getCellSizeX() << " x " 
              << algo.getGrid().getCellSizeY() << "\n";
    std::cout << "  Max occupancy: " << algo.getGrid().getMaxOccupancy() << "\n\n";
    
    algo.run();
    
    std::cout << "\nSimulation completed!\n";
    std::cout << "Total collisions: " << algo.getCollisionCount() << "\n\n";
    
    // Validation
    std::cout << "========== Validation ==========\n";
    if (sys.checkOverlap()) {
        std::cerr << "Final: Overlaps detected!\n";
        return 1;
    }
    std::cout << "Final: No overlaps\n";
    saveConfig(sys, "config_final.dat");
    std::cout << "Saved to config_final.dat\n\n";
    
    // Statistics
    std::cout << "========== Statistics ==========\n";
    std::cout << "  Chains: " << n_chains << "\n";
    std::cout << "  Collisions: " << algo.getCollisionCount() << "\n";
    std::cout << "  Collisions/chain: " 
              << static_cast<double>(algo.getCollisionCount()) / n_chains << "\n\n";
    
    std::cout << "========== TEST PASSED ==========\n";
    std::cout << "Use 'python visualize_config.py' to view results.\n";
    
    return 0;
}
