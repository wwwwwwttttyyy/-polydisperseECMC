#pragma once
#include <iostream>

/**
 * @brief Observer Class - Physical Quantity Measurement and Statistics
 * 
 * Responsibilities:
 *  1. Pressure measurement (X/Y direction separation)
 *     Formula: beta P = rho + (N/V) * sum(deltaIJ) / L_total
 *             = rho * (1 + sum(deltaIJ) / L_total)
 *     Where: deltaIJ = dx_ij - L_event (collision momentum transfer)
 *            L_total = sum(chain_length) (total chain length)
 *            V = box_volume (box volume/area)
 *            rho = N/V (number density)
 *  2. Collision statistics
 *  3. Sampling interval management
 *  4. Future extensions: RDF, energy, configuration output, etc.
 */
class Observer {
public:
    Observer() 
        : collision_count_(0),
          accumulated_delta_x_(0.0),
          accumulated_delta_y_(0.0),
          accumulated_length_x_(0.0),
          accumulated_length_y_(0.0),
          chain_counter_(0),
          configuration_counter_(0) {}
    
    // ========== Update Methods (called by EcmcEngine) ==========
    
    /**
     * @brief Record a collision
     * @param deltaIJ Collision contribution (sqrt(sigma^2 - dy^2) or sqrt(sigma^2 - dx^2))
     * @param direction Direction (0=X, 1=Y)
     */
    void recordCollision(double deltaIJ, int direction) {
        ++collision_count_;
        if (direction == 0) {
            accumulated_delta_x_ += deltaIJ;
        } else {
            accumulated_delta_y_ += deltaIJ;
        }
    }
    
    /**
     * @brief Record event chain length
     * @param chainLength Chain length
     * @param direction Direction (0=X, 1=Y)
     */
    void recordChainLength(double chainLength, int direction) {
        if (direction == 0) {
            accumulated_length_x_ += chainLength;
        } else {
            accumulated_length_y_ += chainLength;
        }
        ++chain_counter_;
    }
    
    /**
     * @brief Reset accumulated quantities for current sampling interval
     */
    void resetSamplingInterval() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        accumulated_length_x_ = 0.0;
        accumulated_length_y_ = 0.0;
        chain_counter_ = 0;
        ++configuration_counter_;
    }
    
    // ========== Query Methods ==========
    
    long long getCollisionCount() const { return collision_count_; }
    long long getChainCounter() const { return chain_counter_; }
    int getConfigurationCounter() const { return configuration_counter_; }
    
    double getAccumulatedDeltaX() const { return accumulated_delta_x_; }
    double getAccumulatedDeltaY() const { return accumulated_delta_y_; }
    double getAccumulatedLengthX() const { return accumulated_length_x_; }
    double getAccumulatedLengthY() const { return accumulated_length_y_; }
    
    /**
     * @brief Calculate pressure factor (1 + sum deltaIJ / L_total)
     * @return Dimensionless pressure factor (divide by volume later for real pressure)
     */
    double calculatePressureFactor() const {
        double total_length = accumulated_length_x_ + accumulated_length_y_;
        double total_delta = accumulated_delta_x_ + accumulated_delta_y_;
        if (total_length > 0) {
            return 1.0 + total_delta / total_length;
        }
        return 1.0;
    }
    
    double calculatePressureFactorX() const {
        if (accumulated_length_x_ > 0) {
            return 1.0 + accumulated_delta_x_ / accumulated_length_x_;
        }
        return 1.0;
    }
    
    double calculatePressureFactorY() const {
        if (accumulated_length_y_ > 0) {
            return 1.0 + accumulated_delta_y_ / accumulated_length_y_;
        }
        return 1.0;
    }
    
    /**
     * @brief Calculate reduced pressure betaP* (normalized with (2r)^2)
     * @param N Number of particles
     * @param volume Box volume
     * @param mean_radius Mean radius
     * @return Reduced pressure
     */
    double calculateReducedPressure(int N, double volume, double mean_radius) const {
        double rho = static_cast<double>(N) / volume;
        double total_length = accumulated_length_x_ + accumulated_length_y_;
        double total_delta = accumulated_delta_x_ + accumulated_delta_y_;
        
        if (total_length <= 0) {
            return rho * 4.0 * mean_radius * mean_radius;  // Ideal gas term
        }
        
        double mean_diameter = 2.0 * mean_radius;
        double normalization = mean_diameter * mean_diameter;
        return (rho + rho * (total_delta / total_length)) * normalization;
    }
    
    /**
     * @brief Calculate X component of reduced pressure
     */
    double calculateReducedPressureX(int N, double volume, double mean_radius) const {
        double rho = static_cast<double>(N) / volume;
        double mean_diameter = 2.0 * mean_radius;
        double normalization = mean_diameter * mean_diameter;
        
        if (accumulated_length_x_ <= 0) {
            return rho * normalization;
        }
        return (rho + rho * (accumulated_delta_x_ / accumulated_length_x_)) * normalization;
    }
    
    /**
     * @brief Calculate Y component of reduced pressure
     */
    double calculateReducedPressureY(int N, double volume, double mean_radius) const {
        double rho = static_cast<double>(N) / volume;
        double mean_diameter = 2.0 * mean_radius;
        double normalization = mean_diameter * mean_diameter;
        
        if (accumulated_length_y_ <= 0) {
            return rho * normalization;
        }
        return (rho + rho * (accumulated_delta_y_ / accumulated_length_y_)) * normalization;
    }
    /**
     * @brief Sample and report pressure
     * @param N Number of particles
     * @param volume Box volume
     * @param mean_radius Mean radius
     * @param sample_interval Sampling interval (number of chains)
     * @return Whether sampling was successful
     */
    bool sampleAndReportPressure(int N, double volume, double mean_radius, long long sample_interval) {
        // Check if enough data available
        if (chain_counter_ < sample_interval) {
            return false;
        }
        
        // Calculate reduced pressure
        double pressure = calculateReducedPressure(N, volume, mean_radius);
        double pressure_x = calculateReducedPressureX(N, volume, mean_radius);
        double pressure_y = calculateReducedPressureY(N, volume, mean_radius);
        
        // Output to terminal
        std::cout << "Pressure [" << configuration_counter_ << "]: "
                  << "betaP* = " << pressure 
                  << " (Px = " << pressure_x << ", Py = " << pressure_y << ")\n";
        
        // TODO: Output to file
        
        // Reset sampling interval
        resetSamplingInterval();
        
        return true;
    }
    
    // ========== Future Extension Interfaces ==========
    // TODO: void recordRDF(const System& sys);
    // TODO: void recordSnapshot(const System& sys, const std::string& filename);
    // TODO: void recordEnergy(double energy);
    
private:
    // Global statistics
    long long collision_count_;      ///< Total collision count (accumulated)
    
    // Sampling interval accumulations (periodically reset)
    double accumulated_delta_x_;     ///< X direction accumulated deltaIJ
    double accumulated_delta_y_;     ///< Y direction accumulated deltaIJ
    double accumulated_length_x_;    ///< X direction accumulated chain length
    double accumulated_length_y_;    ///< Y direction accumulated chain length
    
    // Counters
    long long chain_counter_;        ///< Chain count in current sampling interval
    int configuration_counter_;      ///< Output configuration number
};
