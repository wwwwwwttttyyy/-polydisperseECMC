#define _USE_MATH_DEFINES
#include "initconfig.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
double coslaw(double a, double b, double c) {
    return acos((a * a + b * b - c * c) / (2 * a * b));
}

std::tuple<std::vector<Particle>, double, Box> ini_config(int configtype, int ny, int nx, double target_volume_fraction) {
    std::vector<Particle> data;
    double sigma = 0.0;
    Box box;
    
    if (configtype == 1) {
        double r = 1.0;
        double alp = 0.637556;
        double R = r / alp;
        
        double h1 = sqrt(R * R + 2 * r * R);
        double theta1 = acos(r / (r + R));
        double theta2 = coslaw(r + R, r + R, 2 * R);
        double theta3 = M_PI - (theta1 + theta2);
        double dx = cos(theta3) * (r + R);
        double dy = sin(theta3) * (r + R);
        
        double X = 2 * dx + 2 * r;
        double Y = h1 + R;
        double Lx = X;
        double Ly = Y * 2;
        
        // 计算当前体积分数 - 每个单元有4个小粒子和4个大粒子
        double small_disks_area = 4 * M_PI * r * r;
        double large_disks_area = 4 * M_PI * R * R;
        double total_disks_area = small_disks_area + large_disks_area;
        double cell_area = Lx * Ly;
        double current_eta = total_disks_area / cell_area;
        
        // 调整半径以达到目标体积分数
        if (target_volume_fraction > 0) {
            // 对于 configtype 1，需要保持盒子尺寸不变，只调整粒子半径
            // 计算新的半径来达到目标体积分数
            double target_total_area = target_volume_fraction * Lx * Ly;
            double current_total_area = total_disks_area;
            double area_ratio = target_total_area / current_total_area;
            double radius_scale_factor = sqrt(area_ratio);
            
            r *= radius_scale_factor;
            R *= radius_scale_factor;
            
            // 注意：这里不重新计算几何参数 X, Y, Lx, Ly
            // 保持原来的盒子尺寸，只调整粒子半径
            
            std::cout << "Adjusted radii to achieve target volume fraction: " << target_volume_fraction << std::endl;
            std::cout << "New small radius: " << r << ", large radius: " << R << std::endl;
            std::cout << "Box size remains: " << Lx << " x " << Ly << std::endl;
        }
        
        sigma = 1.0 / (Ly * ny);
        
        if (nx == 0) {
            nx = int(Ly / Lx * ny);
            std::cout << "Automatically set nx=" << nx << " based on ny=" << ny << " and aspect ratio Ly/Lx=" << Ly / Lx << std::endl;
        }
        
        // 计算固定位置（基于初始半径）
        double initial_r = 1.0;
        double initial_R = initial_r / alp;
        
        double x1 = initial_r;
        double y1 = 0;
        double x2 = X / 2 - initial_r;
        double y2 = Y;
        double x3 = X / 2 + initial_r;
        double y3 = Y;
        double x4 = X - initial_r;
        double y4 = 0;
        double x5 = 0;
        double y5 = Y - initial_R;
        double x6 = 0;
        double y6 = initial_R + Y;
        double x7 = X / 2;
        double y7 = initial_R;
        double x8 = X / 2;
        double y8 = 2 * Y - initial_R;
        
        // 使用固定位置，只调整半径
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                data.push_back({x1 + i * X, y1 + j * Ly, r, 1});
                data.push_back({x2 + i * X, y2 + j * Ly, r, 1});
                data.push_back({x3 + i * X, y3 + j * Ly, r, 1});
                data.push_back({x4 + i * X, y4 + j * Ly, r, 1});
                data.push_back({x5 + i * X, y5 + j * Ly, R, 2});
                data.push_back({x6 + i * X, y6 + j * Ly, R, 2});
                data.push_back({x7 + i * X, y7 + j * Ly, R, 2});
                data.push_back({x8 + i * X, y8 + j * Ly, R, 2});
            }
        }
        
        box = {Lx * nx, Ly * ny, Lx, Ly};
        
        // 重新计算最终体积分数 - 每个单元有4个小粒子和4个大粒子
        small_disks_area = 4 * M_PI * r * r;
        large_disks_area = 4 * M_PI * R * R;
        total_disks_area = small_disks_area + large_disks_area;
        cell_area = Lx * Ly;
        double eta = total_disks_area / cell_area;
        
        std::cout << "Initial configuration generated" << std::endl;
        std::cout << "Filling fraction (eta): " << eta << std::endl;
    }
    else if (configtype == 2) {
        // 2D square packing for 1 types, 1 radius, 9 particles in a cell
        double r = 1.0;
        double Lx = 6.0;
        double Ly = 6.0;
        
        // 计算当前体积分数
        double current_eta = 9 * M_PI * r * r / (Lx * Ly);
        
        // 调整半径以达到目标体积分数
        if (target_volume_fraction > 0) {
            double scale_factor = sqrt(target_volume_fraction / current_eta);
            r *= scale_factor;
            std::cout << "Adjusted radius to achieve target volume fraction: " << target_volume_fraction << std::endl;
            std::cout << "New radius: " << r << std::endl;
        }
        
        sigma = 1.0 / (Ly * ny);
        if (nx == 0) {
            nx = static_cast<int>(Ly/Lx * ny);
            std::cout << "Automatically set nx=" << nx << " based on ny=" << ny 
                      << " and aspect ratio Ly/Lx=" << Ly/Lx << std::endl;
        }
        // 定义9个粒子的相对位置
        std::vector<std::pair<double, double>> positions = {
            {1.0, 1.0}, {3.0, 1.0}, {5.0, 1.0},
            {1.0, 3.0}, {3.0, 3.0}, {5.0, 3.0},
            {1.0, 5.0}, {3.0, 5.0}, {5.0, 5.0}
        };
        
        // 添加所有粒子
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (const auto& pos : positions) {
                    data.push_back({pos.first + i * Lx, pos.second + j * Ly, r, 1});
                }
            }
        }
        
        box = {Lx * nx, Ly * ny, Lx, Ly};
        
        // 重新计算最终体积分数
        double eta = 9 * M_PI * r * r / (Lx * Ly);
        std::cout << "Initial configuration generated" << std::endl;
        std::cout << "Filling fraction (eta): " << eta << std::endl;
    }
    else if (configtype == 3) {
        std::vector<std::pair<double, double>> L = {
            {0, 0}, {3.082030, 3.687579}, {4.685883, 6.307616}, {7.767913, 2.620037}
        };
        
        std::vector<std::pair<double, double>> M = {
            {0.412341, 4.588037}, {2.669689, 0.900459}, 
            {5.098223, 1.719578}, {7.355572, 5.407157}
        };
        
        std::vector<std::pair<double, double>> N = {
            {0.919124, 2.363558}, {5.605007, 3.944058}, 
            {2.162906, 6.051137}, {6.848789, 0.256479}
        };
        
        double a = 0.834306042853017;
        double b = 0.651050185882609;
        double R = 1.0 / b;      
        double r = a / b;        
        double smallR = 1.0;     
        double Lx0 = 9.371766;
        double Ly0 = 7.375158;
        
        // 计算当前体积分数
        double total_disks_area_current = 4 * M_PI * R * R + 4 * M_PI * r * r + 4 * M_PI * smallR * smallR;
        double cell_area_current = Lx0 * Ly0;
        double current_eta = total_disks_area_current / cell_area_current;
        
        // 调整半径以达到目标体积分数
        if (target_volume_fraction > 0) {
            double scale_factor = sqrt(target_volume_fraction / current_eta);
            R *= scale_factor;
            r *= scale_factor;
            smallR *= scale_factor;
            
            std::cout << "Adjusted radii to achieve target volume fraction: " << target_volume_fraction << std::endl;
            std::cout << "New radii: R=" << R << ", r=" << r << ", smallR=" << smallR << std::endl;
        }
        
        if (nx == 0) {
            nx = static_cast<int>(Ly0/Lx0 * ny);
            std::cout << "Automatically set nx=" << nx << " based on ny=" << ny 
                      << " and aspect ratio Ly/Lx=" << Ly0/Lx0 << std::endl;
        }
        
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (const auto& pos : L) {
                    data.push_back({pos.first + i * Lx0, pos.second + j * Ly0, R, 3});
                }
                for (const auto& pos : M) {
                    data.push_back({pos.first + i * Lx0, pos.second + j * Ly0, r, 2});
                }
                for (const auto& pos : N) {
                    data.push_back({pos.first + i * Lx0, pos.second + j * Ly0, smallR, 1});
                }
            }
        }
        box = {Lx0 * nx, Ly0 * ny, Lx0, Ly0};
        
        // 重新计算最终体积分数
        double total_disks_area = 0.0;
        for (const auto& p : data) {
            total_disks_area += M_PI * p.radius * p.radius;
        }
        
        double cell_area = box.width * box.height;
        double eta = total_disks_area / cell_area;
        
        sigma = 1.0 / (Ly0 * ny);
        
        std::cout << "Initial configuration generated" << std::endl;
        std::cout << "Filling fraction (eta): " << eta << std::endl;
    }
    return {data, sigma, box};
}

void write_lammps_data(const std::vector<Particle>& data, const Box& box) {
    std::string filename = "2d_2packing.data";
    double max_radius = 0.0;
    
    for (const auto& p : data) {
        if (p.radius > max_radius) max_radius = p.radius;
    }
    
    std::ofstream f(filename);
    
    f << "LAMMPS data file\n\n";
    f << data.size() << " atoms\n";
    f << "2 atom types\n";
    f << "0.0 " << box.width << " xlo xhi\n";
    f << "0.0 " << box.height << " ylo yhi\n";
    f << -1.01 * max_radius << " " << 1.01 * max_radius << " zlo zhi\n\n";
    
    f << "Atoms\n\n";
    int atom_id = 1;
    for (const auto& p : data) {
        f << atom_id << " " << p.type << " " << p.radius * 2 << " " << 1.0 / (M_PI * p.radius * p.radius) << " ";
        f << p.x << " " << p.y << " 0.0\n";
        atom_id++;
    }
    
    std::cout << "2D data file '" << filename << "' generated." << std::endl;
}

void write_ecmc_config(const std::vector<Particle>& data, const Box& box, const std::string& filename) {
    std::ofstream f(filename);
    
    // ECMC格式：第一行是盒子尺寸 Lx Ly
    f << box.width << " " << box.height << "\n";
    
    // 后续行是每个粒子的坐标和半径：x y radius
    for (const auto& p : data) {
        f << p.x << " " << p.y << " " << p.radius << "\n";
    }
    
    std::cout << "ECMC configuration file '" << filename << "' generated." << std::endl;
    std::cout << "Box size: " << box.width << " x " << box.height << std::endl;
    std::cout << "Number of particles: " << data.size() << std::endl;
    
    // 计算并显示体积分数
    double total_area = 0.0;
    for (const auto& p : data) {
        total_area += M_PI * p.radius * p.radius;
    }
    double volume_fraction = total_area / (box.width * box.height);
    std::cout << "Volume fraction (packing fraction): " << volume_fraction << std::endl;
}
