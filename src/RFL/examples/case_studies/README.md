# RFL Case Studies and Specialised Examples

This directory contains advanced use cases, specific physics simulations, and legacy thesis code for the Random Fuzzy Library. 

Basic usage examples are located in `examples/cpp/main.cpp` and `examples/python/main.py`. The code in this directory applies RFL to specialised research problems.

## Contents

- **`type_13_simulation/`**: C++ simulation for the Type (1, 3) Dirac operator using the RFL architecture. Includes eigenvalue recording and Metropolis evolution.
- **`hmc_tuning.cpp`**: Demonstrates parameter tuning (such as `dt` and target acceptance rate) for the Hamiltonian Monte Carlo (HMC) algorithm using the legacy `Geom24` API.
- **`mauro_thesis_mmc.cpp`**: Historical simulation code from Mauro's thesis work, demonstrating the legacy `Geom24` API running a Metropolis Monte Carlo (MMC) evolution.
