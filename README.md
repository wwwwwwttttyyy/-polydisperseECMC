# polydisperseECMC

A high-performance C++17 implementation of the Event-Chain Monte Carlo (ECMC) algorithm for simulating hard disk systems in 2D.
This documentation and the comments in the code were generated with assistance from DeepSeek V3.2.

## Features

- **Efficient spatial partitioning**: Adaptive cell grid with automatic optimization
- **Polydisperse systems**: Support for size-distributed particles  
- **Pressure calculation**: Virial pressure estimation with equilibration control
- **Flexible I/O**: Both text and binary configuration formats
- **Modern C++**: Header-only observers, structured bindings, clean architecture

## Quick Start

```cpp
#include "system.h"
#include "ecmc.h"

System sys;
generateSquareLattice(sys, 16, 16, 0.7);  // 256 particles, φ=0.7

AlgorithmPara algo;
algo.set_n_steps(100000);
algo.set_chain_length(1.0);
algo.set_n_equilibration_pressure(10000);  // Skip first 10k chains

EcmcEngine engine(sys, algo, phys);
engine.run();
```

## Build

Requires C++17 compiler:

```bash
cd tests/build
cmake -G Ninja ..
ninja
./test_ecmc
```

## References

This implementation is based on the Event-Chain Monte Carlo algorithm for hard-sphere systems:

1. E.P. Bernard, W. Krauth, and D.B. Wilson  
   *Event-chain Monte Carlo algorithms for hard-sphere systems*  
   Physical Review E 80, 056704 (2009)  
   doi:10.1103/PhysRevE.80.056704

2. M. Michel, S.C. Kapfer, and W. Krauth  
   *Generalized event-chain Monte Carlo: Constructing rejection-free global-balance algorithms*  
   Journal of Chemical Physics 140, 054116 (2014)  
   doi:10.1063/1.4863991

3. W. Krauth  
   *Statistical Mechanics: Algorithms and Computations*  
   Oxford University Press (2006)  
   ISBN: 978-0198515357
