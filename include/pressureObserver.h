#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

/**
 * @brief 压强观测器 - 专门用于压强测量
 * 
 * 物理公式:
 *   βP* = [ρ + ρ·(Σδ/L)] × (2r̄)²
 * 
 * 其中:
 *   ρ = N/V (数密度)
 *   δᵢⱼ = dx - L_event (X方向) 或 dy - L_event (Y方向)
 *   L = n_chains × chain_length (总链长)
 *   (2r̄)² 为归一化因子，r̄ 为平均半径
 * 
 * 使用方法:
 *   1. initialize(...) - 设置系统参数
 *   2. recordChain(direction) - 每条链开始时调用（仅计数）
 *   3. recordCollision(deltaIJ, direction) - 每次碰撞时调用
 *   4. sampleAndReport() - 达到采样间隔时调用
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
          equilibration_samples_(0),
          phase_chain_count_(0) {}
    
    /**
     * @brief 初始化系统参数（仅调用一次）
     * @param N 粒子数
     * @param volume 盒子体积/面积
     * @param mean_radius 平均半径
     * @param chain_length 单条链的长度
     * @param sample_interval 采样间隔（链数）
     */
    void initialize(int N, double volume, double mean_radius, double chain_length, long long sample_interval) {
        N_ = N;
        volume_ = volume;
        mean_radius_ = mean_radius;
        chain_length_ = chain_length;
        sample_interval_ = sample_interval;
    }
    
    /**
     * @brief 记录碰撞贡献
     * @param deltaIJ 碰撞时刻的接触距离
     * @param direction 0=X方向, 1=Y方向
     */
    void recordCollision(double deltaIJ, int direction) {
        if (direction == 0) {
            accumulated_delta_x_ += deltaIJ;
        } else {
            accumulated_delta_y_ += deltaIJ;
        }
    }
    
    /**
     * @brief 记录一条链开始（仅计数）
     * @param direction 0=X方向, 1=Y方向
     * 
     * 应在链开始、碰撞循环之前调用
     */
    void recordChain(int direction) {
        if (direction == 0) {
            ++chain_count_x_;
        } else {
            ++chain_count_y_;
        }
        ++phase_chain_count_;
    }
    
    /**
     * @brief 计算约化压强 βP*
     */
    double calculateReducedPressure() const {
        double rho = static_cast<double>(N_) / volume_;
        long long total_chain_count = chain_count_x_ + chain_count_y_;
        double total_length = total_chain_count * chain_length_;
        double total_delta = accumulated_delta_x_ + accumulated_delta_y_;
        
        double normalization = 4.0 * mean_radius_ * mean_radius_;  // (2r̄)²
        
        if (total_length <= 0) {
            return rho * normalization;  // 理想气体项
        }
        
        return (rho + rho * (total_delta / total_length)) * normalization;
    }
    
    /**
     * @brief 计算 X 方向约化压强
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
     * @brief 计算 Y 方向约化压强
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
     * @brief 采样并输出压强（达到采样间隔时调用）
     * @return 是否成功采样
     * 
     * 注意: 此方法输出**累积平均值**,不会重置累积量
     * 要获得更小的统计误差,应增加采样间隔或使用 clear() 跳过弛豫
     */
    bool sampleAndReport() {
        if (sample_interval_ <= 0) {
            return false;
        }

        long long interval_chains = chain_count_x_ + chain_count_y_;
        if (interval_chains < sample_interval_) {
            return false;
        }
        
        // 计算累积平均压强(从 clear() 或初始化以来的全部数据)
        double pressure = calculateReducedPressure();
        double pressure_x = calculateReducedPressureX();
        double pressure_y = calculateReducedPressureY();
        
        // 记录到历史数据
        pressure_history_.push_back(pressure);
        pressure_x_history_.push_back(pressure_x);
        pressure_y_history_.push_back(pressure_y);
        chain_count_history_.push_back(phase_chain_count_);
        
        std::cout << "Pressure [" << sample_counter_ << "]: "
                  << "betaP* = " << pressure 
                  << " (Px = " << pressure_x << ", Py = " << pressure_y << ")"
                  << " [phase " << phase_chain_count_
                  << " chains, interval " << interval_chains << " chains]\n";
        
        ++sample_counter_;
        reset();
        return true;
    }
    
    /**
     * @brief 清除所有累积数据（用于删除弛豫阶段）
     * 
     * 在弛豫完成、压强稳定后调用此方法，清除之前的数据，
     * 从干净状态开始采样
     */
    void clear() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        chain_count_x_ = 0;
        chain_count_y_ = 0;
        phase_chain_count_ = 0;
        sample_counter_ = 0;
        
        // 标记弛豫结束位置
        equilibration_samples_ = pressure_history_.size();
    }
    
    /**
     * @brief 重置当前采样区间
     */
    void reset() {
        accumulated_delta_x_ = 0.0;
        accumulated_delta_y_ = 0.0;
        chain_count_x_ = 0;
        chain_count_y_ = 0;
    }
    
    // 查询接口
    long long getTotalChainCount() const { return chain_count_x_ + chain_count_y_; }
    int getSampleCounter() const { return sample_counter_; }
    
    /**
     * @brief 输出压强统计到文件
     * @param filename 输出文件名
     */
    void writePressureData(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Cannot open " << filename << " for writing\n";
            return;
        }
        
        // 计算平衡后数据的统计量
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
        
        // 计算方差(仅平衡后)
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
        
        // 写入统计摘要
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
        
        // 写入所有采样数据
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
    // 系统参数（初始化时设置，运行中不变）
    int N_;                          ///< 粒子数
    double volume_;                  ///< 盒子体积
    double mean_radius_;             ///< 平均半径
    double chain_length_;            ///< 单条链长度
    long long sample_interval_;      ///< 采样间隔
    
    // 累积量（每个采样区间重置）
    double accumulated_delta_x_;     ///< X方向累积δ
    double accumulated_delta_y_;     ///< Y方向累积δ
    long long chain_count_x_;        ///< X方向链计数
    long long chain_count_y_;        ///< Y方向链计数
    int sample_counter_;             ///< 采样次数（输出编号）
    
    // 历史数据（用于统计分析和输出）
    std::vector<double> pressure_history_;       ///< 压强历史
    std::vector<double> pressure_x_history_;     ///< X方向压强历史
    std::vector<double> pressure_y_history_;     ///< Y方向压强历史
    std::vector<long long> chain_count_history_; ///< 累积链数历史
    size_t equilibration_samples_;               ///< 弛豫阶段采样次数
    long long phase_chain_count_;    ///< 当前阶段累计链数
};
