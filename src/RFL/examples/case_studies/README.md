# RFL Case Studies and Specialized Examples

This directory contains advanced use cases, specific physics simulations, and legacy thesis code for the Random Fuzzy Library. 

While the official basic usage examples are located in the parent `examples/` directory (such as `examples/cpp/main.cpp` and `examples/python/main.py`), the code here demonstrates how to apply RFL to solve more complex or specialized research problems.

## Contents

- **`type_13_simulation/`**: An extensive object-oriented C++ simulation setup specifically for the Type (1, 3) Dirac Operator using the modern RFL architecture. Includes eigenvalue recording and Metropolis evolution.
- **`hmc_tuning.cpp`**: A script that demonstrates how to tune parameters (like `dt` and `tgt` acceptance rate) for the Hamiltonian Monte Carlo (HMC) algorithm using the legacy `Geom24` API.
- **`mauro_thesis_mmc.cpp`**: Historical simulation code used as part of Mauro's original thesis work, demonstrating the legacy `Geom24` API executing a standard Metropolis Monte Carlo (MMC) evolution.
