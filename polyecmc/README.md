# ECMC - Event-Chain Monte Carlo for Hard Disks

Efficient Event-Chain Monte Carlo (ECMC) simulation program for 2D hard disk systems, supporting polydisperse systems and pressure measurement.

## Features

- ✅ **Efficient Cell Grid Algorithm** - O(1) collision detection
- ✅ **Periodic Boundary Conditions** - Proper handling of particle boundary crossing
- ✅ **Pressure Measurement** - Accurate pressure calculation based on collision accumulation
- ✅ **Polydisperse Systems** - Support for particle radius distribution
- ✅ **Snapshot Output** - Periodic saving of system configurations
- ✅ **Statistical Analysis** - Automatic separation of equilibration and production periods
- ✅ **Binary I/O** - Binary format saving ~9% storage space

## Building

### Requirements

- CMake 3.15+
- C++11 compiler (GCC, Clang, MSVC)
- Ninja (optional, recommended)

### Compilation Steps

```bash
# Windows (PowerShell)
mkdir build
cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja

# Linux/Mac
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

## Usage

### Preparing Initial Configuration File

The program requires an initial configuration file as input. The project provides an example file `config_initial.dat` (256 particles, φ=0.70).

**Text Format**:
```text
33.8958 33.8958              # Line 1: box dimensions Lx Ly
-15.8887 -15.8887 1.01186    # Lines 2+: x y radius (one line per particle)
-15.8887 -13.7702 0.987337
-15.8887 -11.6517 1.01119
...
```

**Note**: Coordinate system has box center as origin, range [-Lx/2, Lx/2) × [-Ly/2, Ly/2)

### 1. Using Parameter File (Recommended)

Create parameter file `params.txt`:

```text
input_file = config_initial.dat
n_chains = 1000000
sample_interval_pressure = 50000
sample_interval_snapshot = 100000
n_equilibration_pressure = 100000
rng_seed = 42
```

Run simulation:

```bash
./ecmc_main -p params.txt
```

### 2. Pure Command Line Mode

```bash
./ecmc_main -i initial.dat -n 1000000 -s 50000 -e 100000
```

### 3. Mixed Mode (Command Line Overrides Parameter File)

```bash
./ecmc_main -p params.txt -n 2000000 -s 100000
```

## Parameter Description

| Parameter | Short Option | Default | Description |
|-----------|--------------|---------|-------------|
| `input_file` | `-i` | *Required* | Initial configuration file |
| `output_file` | `-o` | `config_final.dat` | Output configuration file |
| `n_chains` | `-n` | 100000 | Total number of event chains |
| `chain_length` | `-l` | 1.0 | Chain length (box size multiples) |
| `sample_interval_pressure` | `-s` | 10000 | Pressure sampling interval |
| `sample_interval_snapshot` | `-S` | 0 | Snapshot output interval (0=disabled) |
| `n_equilibration_pressure` | `-e` | 0 | Number of chains for pressure equilibration |
| `rng_seed` | `-r` | *Random* | Random number seed |
| `binary_output` | `-b` | false | Use binary output |

## Output Files

### 1. Final Configuration (`config_final.dat`)

Contains particle positions and radii at simulation end.

**Text Format**:
```text
256         # Number of particles
33.8958 33.8958    # Box dimensions
-3.142 5.678 1.000  # x y radius
...
```

**Binary Format**: Use `-b` option, saves ~9% space.

### 2. Pressure Data (`pressure_data.txt`)

Contains complete pressure statistics:

```text
# Production Statistics (after equilibration):
#   betaP* = 9.28974 +/- 0.0483076
#   Px     = 9.31525 +/- 0.0632583
#   Py     = 9.26423 +/- 0.0548891
#
# Sample  Chains  betaP*   Px       Py       Phase
0  50000  9.5959  9.6123  9.5795  equilibration
1  100000 9.3621  9.4012  9.3230  production
...
```

### 3. Snapshot Files (`snapshot_*.dat`)

Periodically saved configuration files, filename format: `snapshot_00000100000.dat`

Enabled with `-S` option, e.g., `-S 100000` saves every 100,000 chains.

## Pressure Calculation Formula

This program uses collision-based pressure calculation:

$$
\beta P^* = \left[\rho + \rho \cdot \frac{\sum_i \delta_i}{L}\right] \times (2\bar{r})^2
$$

Where:
- $\rho = N/V$ - Number density
- $\delta_i$ - Free path reduction for the $i$-th collision
- $L = n_{\text{chains}} \times l_{\text{chain}}$ - Total chain length
- $\bar{r}$ - Average particle radius

## Testing

```bash
cd tests/build
cmake .. -G Ninja
ninja
ctest
```

Tests include:
- Single particle movement
- Collision detection
- Periodic boundary conditions
- High-density system stability
- Polydisperse systems
- I/O format validation

## Performance Recommendations

1. **Equilibration Setting**: `n_equilibration_pressure ≈ 0.1 × n_chains`
2. **Sampling Frequency**: `sample_interval_pressure ≈ 0.05 × n_chains`
3. **Snapshot Frequency**: `sample_interval_snapshot ≈ 0.1 × n_chains` (or disable to save I/O)
4. **Packing Fraction**: φ < 0.75 recommended (avoid glass transition)

## Validation Results

For hard disk system with φ = 0.70:
- **Test Result**: βP* = 9.29 ± 0.05
- **Carnahan-Starling**: βP* ≈ 9.28
- **Relative Error**: < 0.2%

## References

1. Bernard et al., *Phys. Rev. E* **80**, 056704 (2009)
2. Michel et al., *J. Chem. Phys.* **140**, 054116 (2014)
3. Kampmann et al., *J. Chem. Phys.* **147**, 214110 (2017)

## License

MIT License

## Authors

ECMC Implementation for Hard Disks (2024)
