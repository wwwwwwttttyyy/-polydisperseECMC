#include "para.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
void AlgorithmPara::set_algorithm_para(
    long long n_steps_,
    double chain_length_,
    long long sample_interval_pressure_) {
  n_steps = n_steps_;
  chain_length = chain_length_;
  sample_interval_pressure = sample_interval_pressure_;
}

void AlgorithmPara::set_n_steps(long long n_steps_) { n_steps = n_steps_; }
void AlgorithmPara::set_chain_length(double chain_length_) { chain_length = chain_length_; }
void AlgorithmPara::set_sample_interval_pressure(long long sample_interval_pressure_) { sample_interval_pressure = sample_interval_pressure_; }
void AlgorithmPara::set_sample_interval_snapshot(long long sample_interval_snapshot_) { sample_interval_snapshot = sample_interval_snapshot_; }
void AlgorithmPara::set_rng_seed(unsigned long rng_seed_) { rng_seed = rng_seed_; }
void AlgorithmPara::set_n_equilibration_pressure(long long n_equilibration_pressure_) { n_equilibration_pressure = n_equilibration_pressure_; }



long long AlgorithmPara::get_n_steps() const { return n_steps; }
double AlgorithmPara::get_chain_length() const { return chain_length; }
long long AlgorithmPara::get_sample_interval_pressure() const { return sample_interval_pressure; }
long long AlgorithmPara::get_sample_interval_snapshot() const { return sample_interval_snapshot; }
unsigned long AlgorithmPara::get_rng_seed() const { return rng_seed; }
long long AlgorithmPara::get_n_equilibration_pressure() const { return n_equilibration_pressure; }

void PhysicalPara::set_n_atoms(int n_atoms_) { n_atoms = n_atoms_; }
void PhysicalPara::set_VolumeFraction(double VolumeFraction_) 
{   VolumeFraction = VolumeFraction_; 
    density = VolumeFraction * 4.0 / M_PI;
}
void PhysicalPara::set_density(double density_) 
{ density = density_; 
  VolumeFraction = density * M_PI / 4.0;
}
void PhysicalPara::set_box(std::array<double, 2> box_) { box = box_; }

int PhysicalPara::get_n_atoms() const { return n_atoms; }
double PhysicalPara::get_VolumeFraction() const { return VolumeFraction; }
double PhysicalPara::get_density() const { return density; }
std::array<double, 2> PhysicalPara::get_box() const { return box; }
