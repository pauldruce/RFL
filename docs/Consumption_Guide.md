# Consuming RFL in Research Projects

This guide explains how to install RFL locally on your machine and consume it across your research projects, Jupyter notebooks, and external simulations.

---

## 1. Python Package Installation

RFL provides precompiled binary wheels on PyPI and supports local source builds via `pybind11` and `scikit-build-core`.

### Option A: PyPI Binary Installation (Easiest / No C++ Setup Required)
Install standalone precompiled wheels from PyPI:
```bash
pip install pyrfl
```
*(For pre-releases, pass `pip install --pre pyrfl`)*

> [!NOTE]
> Precompiled wheels support:
> * **Linux (`x86_64`):** glibc ≥ 2.28 (Ubuntu ≥ 20.04/22.04, RHEL ≥ 8).
> * **macOS (Apple Silicon `arm64`):** macOS ≥ 14.0 (Sonoma / Sequoia).
> * **macOS (Intel `x86_64`):** macOS ≥ 15.0 (Sequoia).
> * **Python:** CPython 3.9 – 3.13.

### Option B: Local Source Installation
For research consumers compiling from a local Git clone:
```bash
pip install src/RFL
```
*(To update after modifying C++ source code: `pip install --force-reinstall --no-deps src/RFL`)*

### Option C: Editable Installation (For Active Development)
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

### Option A: CMake FetchContent (Recommended)
Add RFL directly to your project's `CMakeLists.txt`:
```cmake
include(FetchContent)
FetchContent_Declare(
    RFL
    GIT_REPOSITORY https://github.com/pauldruce/RFL.git
    GIT_TAG        v0.1.0
)
FetchContent_MakeAvailable(RFL)

add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE RFL::core)
```

### Option B: Direct Linking against Compiled Binary
1. Include the relevant headers:
   ```cpp
   #include "DiracOperator.hpp"
   #include "BarrettGlaser/Metropolis.hpp"
   #include "BarrettGlaser/Action.hpp"
   #include "GslRng.hpp"
   ```
2. Link against `librfl.a`, `armadillo`, and `gsl`.
