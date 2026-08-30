# Consuming RFL in Research Projects

This guide explains how to install RFL locally on your machine and consume it across your research projects, Jupyter notebooks, and external simulations.

---

## 1. Local Python Installation (Recommended for Research)

RFL provides native Python bindings (compiled via pybind11 and scikit-build-core with Armadillo integration).

### Option A: Standard Installation (Recommended / Easiest)
For research consumers who just want to use the library without in-tree build artifacts:
```bash
pip install src/RFL
```
*(To update after modifying C++ source code: `pip install --force-reinstall --no-deps src/RFL`)*

### Option B: Editable Installation (For Active Development)
For developers actively iterating on Python bindings and C++ code:
```bash
pip install -e src/RFL
```

### Quick Python Usage Example
You can import and use `rfl` from any script or Jupyter notebook on your machine:

```python
import rfl
import numpy as np

# 1. Initialise a Dirac Operator of signature (p=1, q=3) and matrix dimension 10
dirac = rfl.DiracOperator(p=1, q=3, dim=10)
print(f"Dirac Operator Type: {dirac.get_type()}")
print(f"Total Matrix Dimension: {dirac.get_matrix_dimension()}")

# 2. Run a Metropolis Monte Carlo step
# Parameters: g_2, g_4, scale, num_steps, seed
metropolis = rfl.Metropolis(g_2=-1.0, g_4=1.0, scale=1.0, num_steps=100, seed=42)
acceptance_rate = metropolis.update_dirac(dirac)
print(f"Metropolis Acceptance Rate: {acceptance_rate * 100:.2f}%")

# 3. Analyse Dirac Eigenvalues
eigenvals = dirac.get_eigenvalues()
print(f"Calculated {len(eigenvals)} eigenvalues.")
print(f"Min: {np.min(eigenvals):.4f}, Max: {np.max(eigenvals):.4f}, Mean: {np.mean(eigenvals):.4f}")
```

---

## 2. Building & Testing C++ Core

RFL requires `armadillo`, `gsl`, and `cmake`.

### Build & Run Tests
From the root of the `RFL` repository:
```bash
# 1. Configure and build
cmake -B build src/RFL
cmake --build build --target all -j 4

# 2. Run all unit tests
ctest --test-dir build -j 4 --output-on-failure
```

---

## 3. Consuming RFL in C++ Projects

To link RFL into an external C++ executable:
1. Include the relevant headers:
   ```cpp
   #include "DiracOperator.hpp"
   #include "BarrettGlaser/Metropolis.hpp"
   #include "BarrettGlaser/Action.hpp"
   #include "GslRng.hpp"
   ```
2. Link against `libnew_RFL.a`, `armadillo`, and `gsl`.
