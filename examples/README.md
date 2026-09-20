# Example Applications

This directory contains standalone example applications and research case studies for the Random Fuzzy Library (RFL) in C++ and Python.

---

## Directory Overview

| Directory | Language | Description | Key Targets / Entry Points |
| :--- | :--- | :--- | :--- |
| [`cpp/`](cpp/) | C++17 | Core C++ API example. Initialises a Dirac operator, executes Metropolis sampling, and computes eigenvalues. | `build/examples/cpp/main`<br>`examples/cpp/Makefile` |
| [`python/`](python/) | Python | Python bindings example and interactive Jupyter notebook. Demonstrates NumPy integration and MCMC analysis. | `examples/python/main.py`<br>`examples/python/rfl_playground.ipynb` |
| [`case_studies/`](case_studies/) | C++17 | Specialised research simulations, Markov chain parameter tuning, and historical thesis models. | `mauro_thesis_mmc`<br>`hmc_tuning`<br>`Type13Metropolis` |

---

## 1. Cross-Platform CMake Workflow (Recommended)

CMake builds all C++ examples together with the library:

```bash
# 1. Configure and build from repository root
cmake -B build
cmake --build build -j 4

# 2. Run the basic C++ example
./build/examples/cpp/main

# 3. Run case study simulations
./build/examples/case_studies/mauro_thesis_mmc/mauro_thesis_mmc
./build/examples/case_studies/hmc_tuning/hmc_tuning
./build/examples/case_studies/type_13_simulation/Type13Metropolis
```

---

## 2. Makefile & Direct GCC Workflow (Linux & macOS)

Linux and macOS developers can build examples directly without CMake:

### Option A: Use the Provided Makefile
Navigate to [`examples/cpp/`](cpp/) and use standard `make`:

```bash
cd examples/cpp

# Build and execute
make run

# Clean build artifacts
make clean
```

### Option B: Direct GCC / Clang Command
After compiling `librfl_core.a` in `build/src/core`, compile directly with `g++`:

```bash
g++ -std=c++17 -O3 examples/cpp/main.cpp \
    -Isrc/core \
    -Lbuild/src/core \
    -lrfl_core -larmadillo -lgsl -lgslcblas \
    -o main_gcc

./main_gcc
```

---

## 3. Python & Jupyter Workflow

### Option A: Standard Python CLI
Install the local Python package, then execute the example script:

```bash
pip install .
python3 examples/python/main.py
```

Or run with `uv`:

```bash
uv run python examples/python/main.py
```

### Option B: Interactive Jupyter Notebook
Explore Dirac spectra, Monte Carlo trajectories, and matrix properties interactively:

```bash
jupyter notebook examples/python/rfl_playground.ipynb
```
