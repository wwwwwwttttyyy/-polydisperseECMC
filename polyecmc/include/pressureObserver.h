#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

/**
 * @brief Pressure observer - specialized for pressure measurement
 * 
 * Physical formula:
 *   βP* = [ρ + ρ·(Σδ/L)] × (2r̄)²
 * 
 * Where:
 *   ρ = N/V (number density)
 *   δᵢⱼ = dx - L_event (X direction) or dy - L_event (Y direction)
 *   L = n_chains × chain_length (total chain length)
 *   (2r̄)² is normalization factor, r̄ is mean radius
 * 
 * Usage:
 *   1. initialize(...) - set system parameters
 *   2. recordChain(direction) - call at start of each chain (count only)
 *   3. recordCollision(deltaIJ, direction) - call for each collision
 *   4. sampleAndReport() - call when sampling interval reached
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
          sample_counter_(0),
          equilibration_samples_(0) {}
    
    /**
     * @brief Initialize system parameters (call only once)
     * @param N number of particles
     * @param volume box volume/area
     * @param mean_radius mean radius
     * @param chain_length single chain length
     * @param sample_interval sampling interval (number of chains)
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
     * @param deltaIJ contact distance at collision
     * @param direction 0=X direction, 1=Y direction
     */
    void recordCollision(double deltaIJ, int direction) {
        if (direction == 0) {
            accumulated_delta_x_ += deltaIJ;
        } else {
            accumulated_delta_y_ += deltaIJ;
        }
    }
    
    /**
     * @brief Record start of a chain (count only)
     * @param direction 0=X direction, 1=Y direction
     * 
     * Should be called at chain start, before collision loop
     */
    void recordChain(int direction) {
        if (direction == 0) {
            ++chain_count_x_;
        } else {
            ++chain_count_y_;
        }
    }
    
    /**
     * @brief Calculate reduced pressure βP*
     */
    double calculateReducedPressure() const {
        double rho = static_cast<double>(N_) / volume_;
        long long total_chain_count = chain_count_x_ + chain_count_y_;
        double total_length = total_chain_count * chain_length_;
        double total_delta = accumulated_delta_x_ + accumulated_delta_y_;
        
        double normalization = 4.0 * mean_radius_ * mean_radius_;  // (2r̄)²
        
        if (total_length <= 0) {
            return rho * normalization;  // Ideal gas term
        }
        
        return (rho + rho * (total_delta / total_length)) * normalization;
    }
    
    /**
     * @brief Calculate X direction reduced pressure
     */
    double calculateReducedPressureX() const {
        double rho = static_cast<double>(N_) / volume_;
        double normalization = 4.0 * mean_radius_ * mean_radius_;
        double length_x = chain_count_x_ * chain_length_;
        
        if (length_x <= 0) {
            return rho * normalization;
        }
        
        return (rho + rho * (accumulated_delta_x_ / length_x)) * normalization;
    }
    
    /**
     * @brief Calculate Y direction reduced pressure
     */
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
     * @brief Sample and report pressure (call when sampling interval reached)
     * @return whether sampling was successful
     * 
     * Note: This method outputs **cumulative average**, does not reset accumulations
     * To get smaller statistical error, increase sampling interval or use clear() to skip equilibration
     */
    bool sampleAndReport() {
        long long total_chains = chain_count_x_ + chain_count_y_;
        if (total_chains < sample_interval_) {
            return false;
        }
        
        // Calculate cumulative average pressure (all data since clear() or initialization)
        double pressure = calculateReducedPressure();
        double pressure_x = calculateReducedPressureX();
        double pressure_y = calculateReducedPressureY();
        
        // Record to history data
        pressure_history_.push_back(pressure);
        pressure_x_history_.push_back(pressure_x);
        pressure_y_history_.push_back(pressure_y);
        chain_count_history_.push_back(total_chains);
        
        std::cout << "Pressure [" << sample_counter_ << "]: "
                  << "betaP* = " << pressure 
                  << " (Px = " << pressure_x << ", Py = " << pressure_y << ")"
                  << " [total " << total_chains << " chains]\n";
        
        ++sample_counter_;
        // Note: Do not call reset(), continue accumulating to reduce variance
        return true;
    }
    
    /**
     * @brief Clear all accumulated data (used to remove equilibration phase)
     * 
     * Call this method after equilibration completes and pressure stabilizes,
     * clears previous data, starts sampling from clean state
     */
    void clear() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        chain_count_x_ = 0;
        chain_count_y_ = 0;
        sample_counter_ = 0;
        
        // Mark equilibration end position
        equilibration_samples_ = pressure_history_.size();
    }
    
    /**
     * @brief Reset current sampling interval
     */
    void reset() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        chain_count_x_ = 0;
        chain_count_y_ = 0;
        ++sample_counter_;
    }
    
    // Query interface
    long long getTotalChainCount() const { return chain_count_x_ + chain_count_y_; }
    int getSampleCounter() const { return sample_counter_; }
    
    /**
     * @brief Output pressure statistics to file
     * @param filename output filename
     */
    void writePressureData(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Cannot open " << filename << " for writing\n";
            return;
        }
        
        // Calculate statistics for production data
        size_t n_prod = pressure_history_.size() - equilibration_samples_;
        double mean = 0.0, mean_x = 0.0, mean_y = 0.0;
        
        if (n_prod > 0) {
            for (size_t i = equilibration_samples_; i < pressure_history_.size(); ++i) {
                mean += pressure_history_[i];
                mean_x += pressure_x_history_[i];
                mean_y += pressure_y_history_[i];
            }
            mean /= n_prod;
            mean_x /= n_prod;
            mean_y /= n_prod;
        }
        
        // Calculate variance (production only)
        double var = 0.0, var_x = 0.0, var_y = 0.0;
        if (n_prod > 1) {
            for (size_t i = equilibration_samples_; i < pressure_history_.size(); ++i) {
                double diff = pressure_history_[i] - mean;
                double diff_x = pressure_x_history_[i] - mean_x;
                double diff_y = pressure_y_history_[i] - mean_y;
                var += diff * diff;
                var_x += diff_x * diff_x;
                var_y += diff_y * diff_y;
            }
            var /= (n_prod - 1);
            var_x /= (n_prod - 1);
            var_y /= (n_prod - 1);
        }
        
        // Write statistics summary
        ofs << "# ECMC Pressure Data\n";
        ofs << "# Equilibration samples: " << equilibration_samples_ << "\n";
        ofs << "# Production samples: " << n_prod << "\n";
        ofs << "#\n";
        ofs << "# Production Statistics (after equilibration):\n";
        ofs << "#   betaP* = " << mean << " +/- " << std::sqrt(var) << "\n";
        ofs << "#   Px     = " << mean_x << " +/- " << std::sqrt(var_x) << "\n";
        ofs << "#   Py     = " << mean_y << " +/- " << std::sqrt(var_y) << "\n";
        ofs << "#\n";
        ofs << "# Columns: Sample  Chains  betaP*  Px  Py  Phase\n";
        ofs << "#\n";
        
        // Write all sampling data
        for (size_t i = 0; i < pressure_history_.size(); ++i) {
            ofs << i << "  "
                << chain_count_history_[i] << "  "
                << pressure_history_[i] << "  "
                << pressure_x_history_[i] << "  "
                << pressure_y_history_[i] << "  "
                << (i < equilibration_samples_ ? "equilibration" : "production")
                << "\n";
        }
        
        std::cout << "Pressure data saved to " << filename << "\n";
        std::cout << "  Production: betaP* = " << mean << " +/- " << std::sqrt(var) << "\n";
    }

private:
    // System parameters (set at initialization, constant during run)
    int N_;                          ///< Number of particles
    double volume_;                  ///< Box volume
    double mean_radius_;             ///< Mean radius
    double chain_length_;            ///< Single chain length
    long long sample_interval_;      ///< Sampling interval
    
    // Accumulations (reset each sampling interval)
    double accumulated_delta_x_;     ///< X direction accumulated δ
    double accumulated_delta_y_;     ///< Y direction accumulated δ
    long long chain_count_x_;        ///< X direction chain count
    long long chain_count_y_;        ///< Y direction chain count
    int sample_counter_;             ///< Sample count (output number)
    
    // History data (for statistical analysis and output)
    std::vector<double> pressure_history_;       ///< Pressure history
    std::vector<double> pressure_x_history_;     ///< X direction pressure history
    std::vector<double> pressure_y_history_;     ///< Y direction pressure history
    std::vector<long long> chain_count_history_; ///< Cumulative chain count history
    size_t equilibration_samples_;               ///< Equilibration phase sample count
};
