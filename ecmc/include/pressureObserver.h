#pragma once

#include <iostream>

/**
 * @brief Pressure Observer
 * 
 * Formula: betaP* = [rho + rho*(sum_delta/L)] * (2*r_mean)^2
 * where:
 *   rho = N/V
 *   delta_ij = contact distance at collision
 *   L = total chain length
 */
class PressureObserver {
public:
    PressureObserver() 
        : N_(0),
          volume_(0.0),
          mean_radius_(0.0),
          chain_length_(0.0),
          sample_interval_(0),
          accumulated_delta_x_(0.0),
          accumulated_delta_y_(0.0),
          chain_count_x_(0),
          chain_count_y_(0),
          sample_counter_(0) {}
    
    /**
     * @brief Initialize system parameters (call once)
     */
    void initialize(int N, double volume, double mean_radius, double chain_length, long long sample_interval) {
        N_ = N;
        volume_ = volume;
        mean_radius_ = mean_radius;
        chain_length_ = chain_length;
        sample_interval_ = sample_interval;
    }
    
    /**
     * @brief Record collision contribution
     */
    void recordCollision(double deltaIJ, int direction) {
        if (direction == 0) {
            accumulated_delta_x_ += deltaIJ;
        } else {
            accumulated_delta_y_ += deltaIJ;
        }
    }
    
    /**
     * @brief Record chain start (count only)
     */
    void recordChain(int direction) {
        if (direction == 0) {
            ++chain_count_x_;
        } else {
            ++chain_count_y_;
        }
    }
    
    /**
     * @brief Calculate reduced pressure betaP*
     */
    double calculateReducedPressure() const {
        double rho = static_cast<double>(N_) / volume_;
        long long total_chain_count = chain_count_x_ + chain_count_y_;
        double total_length = total_chain_count * chain_length_;
        double total_delta = accumulated_delta_x_ + accumulated_delta_y_;
        
        double normalization = 4.0 * mean_radius_ * mean_radius_;  // (2*r_mean)^2
        
        if (total_length <= 0) {
            return rho * normalization;  
        }
        
        return (rho + rho * (total_delta / total_length)) * normalization;
    }
    
    double calculateReducedPressureX() const {
        double rho = static_cast<double>(N_) / volume_;
        double normalization = 4.0 * mean_radius_ * mean_radius_;
        double length_x = chain_count_x_ * chain_length_;
        
        if (length_x <= 0) {
            return rho * normalization;
        }
        
        return (rho + rho * (accumulated_delta_x_ / length_x)) * normalization;
    }
    
    double calculateReducedPressureY() const {
        double rho = static_cast<double>(N_) / volume_;
        double normalization = 4.0 * mean_radius_ * mean_radius_;
        double length_y = chain_count_y_ * chain_length_;
        
        if (length_y <= 0) {
            return rho * normalization;
        }
        
        return (rho + rho * (accumulated_delta_y_ / length_y)) * normalization;
    }
    
    /**
     * @brief Sample and report cumulative average pressure
     * @return True if sampling succeeded
     * 
     * Note: Reports cumulative average, does not reset accumulators
     */
    bool sampleAndReport() {
        long long total_chains = chain_count_x_ + chain_count_y_;
        if (total_chains < sample_interval_) {
            return false;
        }
        
        double pressure = calculateReducedPressure();
        double pressure_x = calculateReducedPressureX();
        double pressure_y = calculateReducedPressureY();
        
        std::cout << "Pressure [" << sample_counter_ << "]: "
                  << "betaP* = " << pressure 
                  << " (Px = " << pressure_x << ", Py = " << pressure_y << ")"
                  << " [total " << total_chains << " chains]\n";
        
        ++sample_counter_;
        return true;
    }
    
    /**
     * @brief Clear all accumulated data (e.g., after equilibration)
     */
    void clear() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        chain_count_x_ = 0;
        chain_count_y_ = 0;
        sample_counter_ = 0;
    }
    
    // Reset current sampling interval (unused with cumulative averaging)
    void reset() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        chain_count_x_ = 0;
        chain_count_y_ = 0;
        ++sample_counter_;
    }
    
    long long getTotalChainCount() const { return chain_count_x_ + chain_count_y_; }
    int getSampleCounter() const { return sample_counter_; }

private:
    // System parameters (set once at initialization)
    int N_;                          ///< Particle count
    double volume_;                  ///< Box volume
    double mean_radius_;             ///< Mean radius
    double chain_length_;            ///< Chain length
    long long sample_interval_;      ///< Sampling interval
    
    // Accumulators (cumulative over all chains)
    double accumulated_delta_x_;     ///< X-direction delta sum
    double accumulated_delta_y_;     ///< Y-direction delta sum
    long long chain_count_x_;        ///< X-direction chain count
    long long chain_count_y_;        ///< Y-direction chain count
    int sample_counter_;             ///< Sample counter
};
