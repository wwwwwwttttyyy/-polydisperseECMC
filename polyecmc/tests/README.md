# ECMC Tests

This directory contains a simple demonstration of the Event-Chain Monte Carlo (ECMC) simulation engine.

## Building and Running

```bash
cd tests
mkdir build
cd build
cmake ..
cmake --build .
```

Run the test:
```bash
./test_ecmc          # Linux/macOS
.\test_ecmc.exe      # Windows
```

## What the Test Does

The `test_ecmc.cpp` example demonstrates:

1. **System Initialization**: Creates a 16×16 square lattice of hard disks
2. **ECMC Simulation**: Runs 100,000 event-driven chains
3. **Validation**: Checks for particle overlaps before and after simulation
4. **Output**: Saves initial and final configurations to files

Expected output:
- `config_initial.dat`: Initial square lattice configuration
- `config_final.dat`: Equilibrated liquid-state configuration

## Visualization

To visualize the configurations, use the provided Python script:

```bash
python visualize_config.py
```

This will create a comparison plot showing the initial ordered lattice and the final equilibrated state.

**Requirements**: `matplotlib` and `numpy`
```bash
pip install matplotlib numpy
```

## Configuration Format

Each `.dat` file contains:
- First line: `Lx Ly` (box dimensions)
- Subsequent lines: `x y radius` (particle coordinates and radius)

The coordinate system uses a centered box: `[-L/2, L/2)` with periodic boundary conditions.

## Expected Results

For the default configuration (256 particles, φ=0.6):
- Simulation should complete in ~1 second
- No overlaps detected in initial or final configuration
- Approximately 100-200k total collisions
- Visual transformation from ordered lattice → disordered liquid

## Customization

Edit the configuration section in `test_ecmc.cpp`:

```cpp
const int nx = 16;           // Grid size (nx × ny particles)
const int ny = 16;
const double phi = 0.6;      // Packing fraction
const int n_chains = 100000; // Number of chains
```

Higher packing fractions (φ > 0.7) will have more collisions and take longer to equilibrate.
