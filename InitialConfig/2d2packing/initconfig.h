#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#include <string>

struct Particle {
    double x;
    double y;
    double radius;
    int type;
};

struct Box {
    double width;
    double height;
    double cellWidth;
    double cellHeight;
};

std::tuple<std::vector<Particle>, double, Box> ini_config(int configtype, int ny, int nx = 0, double target_volume_fraction = 0.0);
double coslaw(double a, double b, double c);
void write_lammps_data(const std::vector<Particle>& data, const Box& box);
void write_ecmc_config(const std::vector<Particle>& data, const Box& box, const std::string& filename = "config_initial.dat");
