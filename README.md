# Random Fuzzy Library (RFL)

[![PyPI version](https://img.shields.io/pypi/v/pyrfl.svg?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/pyrfl/)
[![Python versions](https://img.shields.io/pypi/pyversions/pyrfl.svg?logo=python&logoColor=white)](https://pypi.org/project/pyrfl/)
[![GitHub Release](https://img.shields.io/github/v/release/pauldruce/RFL?color=informational&logo=github)](https://github.com/pauldruce/RFL/releases)
[![Build and Test](https://github.com/pauldruce/RFL/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/pauldruce/RFL/actions/workflows/build_and_test.yml)
[![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20WSL2-lightgrey.svg)](https://pypi.org/project/pyrfl/)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg?logo=c%2B%2B)](https://en.cppreference.com/w/cpp/17)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pauldruce/RFL/blob/main/src/RFL/examples/python/rfl_playground.ipynb)

> *"Continuous points are an approximation. Quantise spacetime with finite matrices."*

RFL is a high-performance C++ scientific library with Python bindings for Markov Chain Monte Carlo (MCMC) simulations of **Finite Noncommutative Geometries (Finite NCGs)** and **Random Fuzzy Spaces**.

It simulates the Barrett-Glaser spectral action:

$$S(D) = g_2 \mathrm{Tr}(D^2) + g_4 \mathrm{Tr}(D^4)$$

where $g_2, g_4 \in \mathbb{R}$ and $D$ is a finite Dirac operator. See §[Background & Citations](#background--citations) below for foundational academic papers.

---

## Quick Start (Python)

Precompiled binary wheels vendor dynamic linear algebra dependencies (`OpenBLAS`, `Armadillo`, `GSL`). No C++ compiler or local dependencies are required.

### 1. Installation

Install the package from PyPI:

```bash
pip install pyrfl
```

> [!NOTE]
> The PyPI distribution package name is **`pyrfl`**. The Python module name is **`rfl`** (`import rfl`).

### 2. Run a Simulation

Run this 30-second simulation:

```python
import numpy as np
import rfl

# Initialise a Dirac operator with signature (p=1, q=3) and matrix dimension 10.
dirac = rfl.DiracOperator(p=1, q=3, dim=10)

# Configure and run the Metropolis algorithm.
# Parameters: g_2, g_4, scale, num_steps, seed
metropolis = rfl.Metropolis(g_2=-1.0, g_4=1.0, scale=1.0, num_steps=100, seed=42)
acceptance_rate = metropolis.update_dirac(dirac)
print(f"Acceptance Rate: {acceptance_rate * 100:.2f}%")

# Extract eigenvalues directly as a NumPy array (zero-copy memory mapping).
eigenvals = dirac.get_eigenvalues()
print(f"Calculated {len(eigenvals)} eigenvalues.")
print(f"Spectrum range: [{np.min(eigenvals):.2f}, {np.max(eigenvals):.2f}]")
```

### 3. Interactive Playground

Explore simulations interactively in the Jupyter notebook playground:
* Run locally: [`src/RFL/examples/python/rfl_playground.ipynb`](src/RFL/examples/python/rfl_playground.ipynb)
* Run in browser: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pauldruce/RFL/blob/main/src/RFL/examples/python/rfl_playground.ipynb)

---

## Quick Start (C++)

Add RFL directly to your `CMakeLists.txt` file using `FetchContent`:

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

Include the headers in your C++ code:

```cpp
#include "DiracOperator.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "BarrettGlaser/Action.hpp"
#include "GslRng.hpp"
```

---

## Supported Platforms & Delivery

RFL supports two delivery methods:

### 1. Precompiled Python Wheels (Binary Delivery)

Precompiled wheels vendor dynamic linear algebra dependencies (`openblas`, `gsl`, `armadillo`) directly inside the package:

| Platform | Architecture | Package Standard | Minimum OS Baseline | Python Coverage |
| :--- | :--- | :--- | :--- | :--- |
| **Linux** | `x86_64` | `manylinux_2_28` | glibc ≥ 2.28 (Ubuntu ≥ 20.04/22.04, RHEL ≥ 8) | 3.9 – 3.13 |
| **macOS (Apple Silicon)** | `arm64` | `macosx_14_0_arm64` | macOS ≥ 14.0 (Sonoma / Sequoia) | 3.9 – 3.13 |
| **macOS (Intel)** | `x86_64` | `macosx_15_0_x86_64` | macOS ≥ 15.0 (Sequoia) | 3.9 – 3.13 |
| **Windows** | `x86_64` | WSL2 | Ubuntu 22.04 on WSL2 (Native MSVC in v0.3.0) | 3.9 – 3.13 |

### 2. C++ Source Builds & CMake `FetchContent` (Source Delivery)

When compiling RFL from source or linking via CMake `FetchContent`, RFL compiles against local toolchains:

| Dependency | Minimum Version | Notes |
| :--- | :--- | :--- |
| **C++ Standard** | **C++17** | `std::optional`, `std::variant`, structured bindings |
| **Compilers** | **GCC ≥ 10, Clang ≥ 11, Apple Clang ≥ 13, MSVC ≥ 2019** | Conforming C++17 compilers |
| **CMake** | **≥ 3.20** | Modern target export and `FetchContent` syntax |
| **Armadillo** | **≥ 11.4.0** | High-performance matrix mathematics |
| **GSL** | **≥ 2.6** | GNU Scientific Library random number generators |
| **BLAS / LAPACK** | **OpenBLAS / Accelerate / MKL** | Linear algebra backend for Armadillo |

---

## Building from Source (C++)

To contribute to RFL or build the C++ library locally:

### 1. Install Dependencies

* **macOS:** `brew install cmake gsl armadillo openblas`
* **Ubuntu / Debian:** `sudo apt-get install cmake libgsl-dev libarmadillo-dev libopenblas-dev`
* **Other platforms:** Follow installation instructions provided by CMake, GSL, and Armadillo.

### 2. Build and Test

```bash
git clone https://github.com/pauldruce/RFL.git
cd RFL
cmake -B build src/RFL
cmake --build build --target all -j 4
ctest --test-dir build -j 4 --output-on-failure
```

### 3. Available CMake Targets

* `rfl_core` (aliases: `RFL::core`, `RFL::rfl`): The modern C++17 library.
* `rfl_legacy` (alias: `RFL::legacy`): The preserved historical codebase for baseline comparison.
* `playground`: Standalone C++ experimentation target.

To inspect all available build targets:

```bash
cmake --build build --target help
```

---

## Development & Feature Request Workflow

RFL follows a **3-tier research-driven development workflow** to keep friction near zero while maintaining architectural integrity:

```
1. GitHub Issues (30s)    ──>  Catchment basin for quick bugs, missing functions, and friction notes.
2. EPs in docs/eps/       ──>  Enhancement Proposals for major architectural milestones & new physics.
3. Milestones & PRs       ──>  Execution trackers and deliverable PR slices.
```

### 1. Requesting Small Features or Reporting Friction (GitHub Issues)
When running research experiments (in Jupyter notebooks or scripts) and hitting a missing feature or bug:
* Open a quick **GitHub Issue** (e.g. `Clifford gamma matrices should be accessible directly in Python`).
* Tag with `enhancement`, `bug`, or `research-need`.

### 2. Major Architecture & Physics Upgrades (Enhancement Proposals - EPs)
When several related issues point to a major subsystem upgrade (such as Value Semantics or Fermion Pfaffian Actions):
* Author an **Enhancement Proposal** in `docs/eps/` (e.g. [EP-1](docs/eps/ep-1-core-architecture-modernisation.md)).
* EPs capture **Research Scenarios**, **Requirements & Invariants Table**, and **Architecture Decision Records (ADRs)** with explicit trade-offs.
* EPs are committed in git alongside the code for permanent versioned traceability.

### 3. Milestone Scheduling & PR Delivery
* When an EP is approved, assign the issues to a **GitHub Milestone** (e.g. `v0.3.0: Core Modernisation (EP-1)`).
* The EP delivery plan is converted into discrete GitHub Issues assigned to the milestone.
* PRs reference their corresponding issue (`Closes #12`), enabling automatic milestone progress tracking and issue closure upon merge.

---

## Documentation

Documentation for the software architecture, release lifecycle, and controlled vocabulary is available in the `docs/` directory:
- [Target Architecture Guide](docs/Architecture.md)
- [Release Process & Notes Guide](docs/Release_Process.md)
- [Controlled Vocabulary & Glossary](docs/Glossary.md)
- [Consuming RFL in Research](docs/Consumption_Guide.md)
- [Enhancement Proposals](docs/eps/)
- [v0.1.0 Release Notes](docs/releases/v0.1.0.md)

---

## Background & Citations

Random Finite NCGs are an active area of academic research in Quantum Gravity and Mathematical Physics. Foundational academic papers include:

1. John W. Barrett and Lisa Glaser. (2016). *Monte Carlo simulations of random non-commutative geometries*. Journal of Physics A: Mathematical and Theoretical 49, 24: 245001. [https://doi.org/10.1088/1751-8113/49/24/245001](https://doi.org/10.1088/1751-8113/49/24/245001)
2. Lisa Glaser. (2017). *Scaling behaviour in random non-commutative geometries*. Journal of Physics A: Mathematical and Theoretical 50, 27: 275201. [https://doi.org/10.1088/1751-8121/aa7424](https://doi.org/10.1088/1751-8121/aa7424)
3. John W. Barrett, Paul Druce, and Lisa Glaser. (2019). *Spectral estimators for finite non-commutative geometries*. Journal of Physics A: Mathematical and Theoretical 52, 27: 275203. [https://doi.org/10.1088/1751-8121/ab22f8](https://doi.org/10.1088/1751-8121/ab22f8)
4. Lisa Glaser and Abel Stern. (2020). *Understanding truncated non-commutative geometries through computer simulations*. Journal of Mathematical Physics 61, 3: 033507. [https://doi.org/10.1063/1.5131864](https://doi.org/10.1063/1.5131864)
5. Lisa Glaser and Abel B. Stern. (2021). *Reconstructing manifolds from truncations of spectral triples*. Journal of Geometry and Physics. [https://doi.org/10.1016/j.geomphys.2020.103921](https://doi.org/10.1016/j.geomphys.2020.103921)
6. John W. Barrett. (2015). *Matrix geometries and fuzzy spaces as finite spectral triples*. Journal of Mathematical Physics 56, 082301. [https://doi.org/10.1063/1.4927224](https://doi.org/10.1063/1.4927224)

### BibTeX

```bibtex
@article{Barrett_2016,
  author = {Barrett, John W. and Glaser, Lisa},
  title = {Monte Carlo simulations of random non-commutative geometries},
  journal = {Journal of Physics A: Mathematical and Theoretical},
  volume = {49},
  number = {24},
  pages = {245001},
  year = {2016},
  doi = {10.1088/1751-8113/49/24/245001}
}

@article{Barrett_2019,
  author = {Barrett, John W. and Druce, Paul and Glaser, Lisa},
  title = {Spectral estimators for finite non-commutative geometries},
  journal = {Journal of Physics A: Mathematical and Theoretical},
  volume = {52},
  number = {27},
  pages = {275203},
  year = {2019},
  doi = {10.1088/1751-8121/ab22f8}
}
```