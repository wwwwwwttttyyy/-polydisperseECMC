#include "system.h"
#include "ecmc.h"
#include "para.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <algorithm>

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

// Check if the entire system has no overlaps
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
                std::cout << "  Radius i: " << sys.radius[i] << ", Radius j: " << sys.radius[j] << std::endl;
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

// Sample radii from inverse power law distribution
// P(r) ∝ r^(-α), using α = 3 here
std::vector<double> sampleInversePowerLaw(int N, double rmin, double rmax, int seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    // Inverse power law distribution: P(r) ∝ r^(-α)
    // CDF: F(r) = (r^(1-α) - rmin^(1-α)) / (rmax^(1-α) - rmin^(1-α))
    // Inverse sampling: r = [u * (rmax^(1-α) - rmin^(1-α)) + rmin^(1-α)]^(1/(1-α))
    
    double alpha = 3.0;  // Power exponent
    double power = 1.0 - alpha;  // = -2
    
    double rmin_pow = std::pow(rmin, power);
    double rmax_pow = std::pow(rmax, power);
    double denom = rmax_pow - rmin_pow;
    
    std::vector<double> radii;
    radii.reserve(N);
    
    for (int i = 0; i < N; ++i) {
        double u = uniform(rng);
        double r_pow = u * denom + rmin_pow;
        double r = std::pow(r_pow, 1.0 / power);
        radii.push_back(r);
    }
    
    return radii;
}

// Generate polydisperse square lattice initial configuration
void generatePolydisperseSquareLattice(System& sys, int nx, int ny, 
                                       double rmin, double rmax, 
                                       double packing_fraction) {
    sys.clear();
    
    int N = nx * ny;
    
    // Generate radius distribution
    auto radii = sampleInversePowerLaw(N, rmin, rmax);
    
    // Calculate mean radius
    double r_mean = 0.0;
    for (double r : radii) {
        r_mean += r;
    }
    r_mean /= N;
    
    // Calculate lattice constant based on packing fraction
    // packing_fraction = N * <π*r²> / (Lx * Ly)
    // For square lattice, Lx = nx * a, Ly = ny * a
    // Therefore packing_fraction = π * <r²> / a²
    
    double r_sq_mean = 0.0;
    for (double r : radii) {
        r_sq_mean += r * r;
    }
    r_sq_mean /= N;
    
    double a = std::sqrt(M_PI * r_sq_mean / packing_fraction);
    
    // Ensure lattice constant is large enough to avoid overlap of largest particles
    double a_min = 2.0 * rmax * 1.01;  // Leave 1% gap
    if (a < a_min) {
        std::cout << "Warning: Lattice constant too small, adjusting..." << std::endl;
        std::cout << "  Calculated a = " << a << ", minimum required = " << a_min << std::endl;
        a = a_min;
    }
    
    double Lx = nx * a;
    double Ly = ny * a;
    
    sys.boxsize = {Lx, Ly};
    
    std::cout << "Polydisperse lattice parameters:" << std::endl;
    std::cout << "  Radius range: [" << rmin << ", " << rmax << "]" << std::endl;
    std::cout << "  Mean radius: " << r_mean << std::endl;
    std::cout << "  Lattice constant a = " << a << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << std::endl;
    std::cout << "  Box size: [" << Lx << ", " << Ly << "]" << std::endl;
    
    // Generate square lattice positions, randomly assign radii
    std::mt19937 rng(123);  // For position randomization
    std::shuffle(radii.begin(), radii.end(), rng);
    
    int idx = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = (i + 0.5) * a - Lx / 2.0;
            double y = (j + 0.5) * a - Ly / 2.0;
            sys.addParticle(x, y, radii[idx], 0);
            idx++;
        }
    }
    
    // Verify packing fraction
    double area = Lx * Ly;
    double particle_area = 0.0;
    for (size_t i = 0; i < sys.size(); ++i) {
        particle_area += M_PI * sys.radius[i] * sys.radius[i];
    }
    double actual_packing = particle_area / area;
    
    std::cout << "  N = " << sys.size() << " particles" << std::endl;
    std::cout << "  Target packing fraction = " << packing_fraction << std::endl;
    std::cout << "  Actual packing fraction = " << actual_packing << std::endl;
    std::cout << "  Polydispersity = " << (rmax - rmin) / r_mean << std::endl;
}

int main() {
    std::cout << "========== Polydisperse ECMC Test ==========" << std::endl;
    
    // Parameter settings - pressure test
    double rmin = 0.8;
    double rmax = 1.1;
    double packing_fraction = 0.7;
    int nx = 32;  // 32x32 = 1024 particles
    int ny = 32;
    
    std::cout << "System parameters:" << std::endl;
    std::cout << "  Radius range: [" << rmin << ", " << rmax << "]" << std::endl;
    std::cout << "  Target packing fraction = " << packing_fraction << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << " = " << (nx*ny) << " particles" << std::endl;
    
    // Generate polydisperse square lattice initial configuration
    std::cout << "\nGenerating polydisperse square lattice..." << std::endl;
    System sys;
    generatePolydisperseSquareLattice(sys, nx, ny, rmin, rmax, packing_fraction);
    
    int N = sys.size();
    double Lx = sys.boxsize[0];
    double Ly = sys.boxsize[1];
    
    // Check initial configuration
    if (!checkAllNoOverlap(sys)) {
        std::cout << " Initial configuration has overlaps!" << std::endl;
        return 1;
    }
    std::cout << "\n Initial configuration: No overlaps" << std::endl;
    
    // Save initial configuration
    saveConfiguration(sys, "config_initial.dat");
    
    // Configure ECMC parameters
    AlgorithmPara algo;
    algo.set_n_steps(5000000);       // Run 5 million chains (pressure test)
    algo.set_chain_length(1.0);      // Chain length 1.0
    algo.set_sample_interval_pressure(100);
    algo.set_sample_interval_snapshot(1000);
    algo.set_rng_seed(12345);
    
    PhysicalPara phys;
    phys.set_n_atoms(N);
    phys.set_box({Lx, Ly});
    
    // Run ECMC
    std::cout << "\n========== Running ECMC ==========" << std::endl;
    EcmcEngine engine(sys, algo, phys);
    
    engine.run();
    
    // Check final configuration
    std::cout << "\n========== Validation ==========" << std::endl;
    if (!checkAllNoOverlap(sys)) {
        std::cout << "FAILED: Final configuration has overlaps!" << std::endl;
        saveConfiguration(sys, "config_final_FAILED.dat");
        return 1;
    }
    
    std::cout << " Final configuration: No overlaps" << std::endl;
    
    // Save final configuration
    saveConfiguration(sys, "config_final.dat");
    
    // Statistics of radius distribution
    std::cout << "\n========== Radius Statistics ==========" << std::endl;
    double r_min_actual = *std::min_element(sys.radius.begin(), sys.radius.end());
    double r_max_actual = *std::max_element(sys.radius.begin(), sys.radius.end());
    double r_mean = sys.radiusMean;
    
    std::cout << "  Min radius: " << r_min_actual << std::endl;
    std::cout << "  Max radius: " << r_max_actual << std::endl;
    std::cout << "  Mean radius: " << r_mean << std::endl;
    
    std::cout << "\n========== TEST PASSED ==========" << std::endl;
    std::cout << "Initial and final configurations saved." << std::endl;
    std::cout << "Use Python script to visualize." << std::endl;
    
    return 0;
}
