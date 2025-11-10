#pragma once
#include <string>
#include <random>
#include <array>      

class AlgorithmPara {
  public:
  //setting algorithm parameters
  void set_algorithm_para(
      long long n_steps_,
      double chain_length_,
      long long sample_interval_pressure_);

  // Individual setters
  void set_n_steps(long long n_steps_);
  void set_chain_length(double chain_length_);
  void set_sample_interval_pressure(long long sample_interval_pressure_);
  void set_sample_interval_snapshot(long long sample_interval_snapshot_);
  void set_rng_seed(unsigned long rng_seed_);
  void set_n_equilibration_pressure(long long n_equilibration_pressure_);  // Pressure equilibration steps
  
  // Individual getters
  long long get_n_steps() const;
  double get_chain_length() const;
  long long get_sample_interval_pressure() const;
  long long get_sample_interval_snapshot() const;
  unsigned long get_rng_seed() const;
  long long get_n_equilibration_pressure() const;
  private:
    long long  n_steps;
    double chain_length;
    long long sample_interval_pressure;
    long long sample_interval_snapshot;
    unsigned long rng_seed;
    long long n_equilibration_pressure;  // Pressure equilibration steps (default 0 means no skipping)
};

class PhysicalPara {

  public:

  // Individual setters
  void set_n_atoms(int n_atoms_);
  void set_VolumeFraction(double VolumeFraction_);
  void set_density(double density_);
  void set_box(std::array<double, 2> box_);

  // Individual getters
  int get_n_atoms() const;
  double get_VolumeFraction() const;
  double get_density() const;
  std::array<double, 2> get_box() const;

  private:
    int n_atoms;
    double VolumeFraction;
    double density;
    std::array<double, 2> box;
};
