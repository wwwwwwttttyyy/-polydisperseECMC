#include "system.h"
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <iostream>

// ========== Constructors ==========
System::System(size_t n, const std::array<double, 2>& box_) : boxsize(box_) {
    reserve(n);
}


void System::reserve(size_t n) {
    x.reserve(n);
    y.reserve(n);
    radius.reserve(n);
    typeID.reserve(n);
}

void System::clear() {
    x.clear();
    y.clear();
    radius.clear();
    typeID.clear();
    radiusMax = 0.0;
    radiusMean = 0.0;
    radiusMin = std::numeric_limits<double>::infinity();
}


void System::addParticle(double x_, double y_, double r_, int tid) {
    if (r_ <= 0.0) {
        throw std::invalid_argument("Particle radius must be positive");
    }
    
    x.push_back(x_);
    y.push_back(y_);
    radius.push_back(r_);
    typeID.push_back(tid);
    
    // Dynamically update radius statistics
    if (r_ > radiusMax) radiusMax = r_;
    if (r_ < radiusMin) radiusMin = r_;
    
    // Update mean value (incremental update to avoid re-traversal)
    size_t n = radius.size();
    radiusMean = ((n - 1) * radiusMean + r_) / n;
}
void System::removeParticle(size_t index) {
    if (index >= x.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    x.erase(x.begin() + index);
    y.erase(y.begin() + index);
    radius.erase(radius.begin() + index);
    typeID.erase(typeID.begin() + index);
}   



void System::updateMaxMeanRadius() {
    if (radius.empty()) {
        radiusMax = 0.0;
        radiusMean = 0.0;
        radiusMin = std::numeric_limits<double>::infinity();
        return;
    }
    
    // Recalculate all radius statistics (used after batch loading)
    radiusMax = *std::max_element(radius.begin(), radius.end());
    radiusMin = *std::min_element(radius.begin(), radius.end());
    
    double sum = std::accumulate(radius.begin(), radius.end(), 0.0);
    radiusMean = sum / radius.size();
}






bool System::checkOverlap() const {
    size_t N = x.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            double dist = distance(i, j);
            if (dist < (radius[i] + radius[j])) {
                return true; 
            }
        }
    }
    return false;
}

void System::applyPBC(size_t i) {
    if (i >= x.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    x[i] -= boxsize[0] * std::round(x[i] / boxsize[0]);
    y[i] -= boxsize[1] * std::round(y[i] / boxsize[1]);
}

std::array<double, 2> System::displacement(size_t i, size_t j) const {
    if (i >= x.size() || j >= x.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    
    double dx = x[j] - x[i];
    double dy = y[j] - y[i];
    dx -= boxsize[0] * std::round(dx / boxsize[0]);
    dy -= boxsize[1] * std::round(dy / boxsize[1]);
    
    return {dx, dy};
}

double System::distance(size_t i, size_t j) const {
    auto d = displacement(i, j);
    return std::sqrt(d[0]*d[0] + d[1]*d[1]);
}

// ========== File I/O ==========
void System::saveConfig(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Line 1: box dimensions
    ofs << boxsize[0] << " " << boxsize[1] << "\n";
    
    // Subsequent lines: particle data
    for (size_t i = 0; i < x.size(); ++i) {
        ofs << x[i] << " " << y[i] << " " << radius[i] << "\n";
    }
    
    ofs.close();
}

void System::loadConfig(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }
    
    // Clear current data
    clear();
    
    // Read box dimensions
    if (!(ifs >> boxsize[0] >> boxsize[1])) {
        throw std::runtime_error("Failed to read box size from: " + filename);
    }
    
    // Read particle data
    double x_val, y_val, r_val;
    while (ifs >> x_val >> y_val >> r_val) {
        // Apply periodic boundary conditions, ensure coordinates in [-L/2, L/2) range
        x_val -= boxsize[0] * std::round(x_val / boxsize[0]);
        y_val -= boxsize[1] * std::round(y_val / boxsize[1]);
        addParticle(x_val, y_val, r_val);
    }
    
    ifs.close();
    
    // Update radius statistics
    updateMaxMeanRadius();
    
    std::cout << "Loaded " << size() << " particles from " << filename 
              << " (PBC applied, box centered at origin)\n";
}

void System::saveConfigBinary(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open binary file for writing: " + filename);
    }
    
    // Write box dimensions
    ofs.write(reinterpret_cast<const char*>(&boxsize[0]), sizeof(double));
    ofs.write(reinterpret_cast<const char*>(&boxsize[1]), sizeof(double));
    
    // Write particle count
    double N_double = static_cast<double>(x.size());
    ofs.write(reinterpret_cast<const char*>(&N_double), sizeof(double));
    
    // Write particle data (continuous write of x, y, r)
    for (size_t i = 0; i < x.size(); ++i) {
        ofs.write(reinterpret_cast<const char*>(&x[i]), sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&y[i]), sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&radius[i]), sizeof(double));
    }
    
    ofs.close();
}

void System::loadConfigBinary(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Cannot open binary file for reading: " + filename);
    }
    
    // Clear current data
    clear();
    
    // Read box dimensions
    if (!ifs.read(reinterpret_cast<char*>(&boxsize[0]), sizeof(double)) ||
        !ifs.read(reinterpret_cast<char*>(&boxsize[1]), sizeof(double))) {
        throw std::runtime_error("Failed to read box size from binary file: " + filename);
    }
    
    // Read particle count
    double N_double;
    if (!ifs.read(reinterpret_cast<char*>(&N_double), sizeof(double))) {
        throw std::runtime_error("Failed to read particle count from binary file: " + filename);
    }
    size_t N = static_cast<size_t>(N_double);
    
    // Pre-allocate space
    reserve(N);
    
    // Read particle data
    for (size_t i = 0; i < N; ++i) {
        double x_val, y_val, r_val;
        if (!ifs.read(reinterpret_cast<char*>(&x_val), sizeof(double)) ||
            !ifs.read(reinterpret_cast<char*>(&y_val), sizeof(double)) ||
            !ifs.read(reinterpret_cast<char*>(&r_val), sizeof(double))) {
            throw std::runtime_error("Failed to read particle data from binary file: " + filename);
        }
        // Apply periodic boundary conditions, ensure coordinates in [-L/2, L/2) range
        x_val -= boxsize[0] * std::round(x_val / boxsize[0]);
        y_val -= boxsize[1] * std::round(y_val / boxsize[1]);
        addParticle(x_val, y_val, r_val);
    }
    
    ifs.close();
    
    // Update radius statistics
    updateMaxMeanRadius();
    
    std::cout << "Loaded " << size() << " particles from binary file " << filename 
              << " (PBC applied, box centered at origin)\n";
}

// ========== Statistical Calculations ==========
double System::getPackingFraction() const {
    if (x.empty()) {
        return 0.0;
    }
    
    double area = boxsize[0] * boxsize[1];
    if (area <= 0.0) {
        throw std::runtime_error("Invalid box size for packing fraction calculation");
    }
    
    // Calculate total area of all particles
    double total_area = 0.0;
    for (size_t i = 0; i < radius.size(); ++i) {
        total_area += 3.14159265358979323846 * radius[i] * radius[i];
    }
    
    return total_area / area;
}
