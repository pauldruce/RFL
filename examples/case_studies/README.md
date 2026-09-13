# RFL Case Studies and Specialised Examples

This directory contains specialised research simulations, algorithm tuning scripts, and thesis models for the Random Fuzzy Library (RFL).

For basic library usage, see the parent [`examples/`](../) directory.

---

## Contents

### 1. Mauro Thesis Metropolis Simulation (`mauro_thesis_mmc.cpp`)
Historical thesis simulation reproducing Barrett-Glaser Monte Carlo Markov Chain updates on $(2, 0)$ geometries with matrix dimension $N=32$ and coupling $g_2=-3.0$.

* **Build Target:** `mauro_thesis_mmc`
* **Run Command:**
  ```bash
  ./build/examples/case_studies/mauro_thesis_mmc
  ```

---

### 2. Hamiltonian Monte Carlo Parameter Tuning (`hmc_tuning.cpp`)
Demonstrates No-U-Turn / dual-averaging step-size ($\mathrm{d}t$) adaptation targeting an $80\%$ acceptance rate, followed by thermalisation and momentum updates.

* **Build Target:** `hmc_tuning`
* **Run Command:**
  ```bash
  ./build/examples/case_studies/hmc_tuning
  ```
* **Output:** Generates `example_S.txt` (action trajectory) and `example_HL.txt` (energy diagnostics).

---

### 3. Type (1, 3) Dirac Operator Simulation (`type_13_simulation/`)
Simulates a type $(1, 3)$ Dirac operator with matrix dimension $N=10$, recording eigenvalue spectra to structured HDF5 files during Metropolis sweeps.

* **Prerequisites:** Requires the HDF5 C++ library (`brew install hdf5` on macOS or `sudo apt-get install libhdf5-dev` on Debian/Ubuntu).
* **Build Target:** `Type13Metropolis`
* **Run Command:**
  ```bash
  ./build/examples/case_studies/type_13_simulation/Type13Metropolis
  ```
