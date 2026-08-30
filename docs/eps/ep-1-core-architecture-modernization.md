# EP-1: Research-Driven Architecture Modernization for RFL

* **Title:** Research-Driven Architecture Modernization for RFL
* **Author:** Paul Druce
* **Status:** In Discussion
* **Target Version:** RFL v2.0
* **Date:** 2026-08-30

---

## 1. Motivation, Goals & Non-Goals

### 1.1 Problem Statement & Research Context
`RFL` was created to study **Finite Noncommutative Geometries (Random Fuzzy Spaces)** and **Spectral Triples** via Monte Carlo simulations. The core mathematical framework is governed by the Barrett-Glaser action (Barrett & Glaser 2016):
$$S(D) = g_2 \text{Tr}(D^2) + g_4 \text{Tr}(D^4)$$
where the Dirac operator is decomposed into Hermitian and anti-Hermitian matrices ($H_i, L_j$) coupled to Clifford gamma matrices:
$$D = \sum_{i=1}^p \gamma^i \otimes H_i + \sum_{j=1}^q \gamma^{p+j} \otimes L_j$$

As the research program expands into **Fermion Functional Integrals (Barrett 2024)**, **Product Geometries $S_F^2 \otimes \mathcal{F}$ (Barrett 2026)**, and **Random Matrix Spectral Statistics**, the legacy codebase exhibits structural bottlenecks:
1. **Opaque Execution Control:** The legacy `Simulation` class executes a monolithic, blocking loop. Researchers cannot inspect eigenvalue trajectories, stream observables to disk, or apply interactive stopping conditions in Jupyter notebooks.
2. **Entangled Physics & Sampling:** The analytical trace variation formulas $\Delta \text{Tr}(D^2)$ and $\Delta \text{Tr}(D^4)$ are embedded directly inside `Metropolis.cpp`. This prevents reusing the same physics for Hybrid Monte Carlo (HMC) or new multi-term action potentials ($S_6, \text{Pfaffian}$).
3. **Leaky Interface & Const Violation:** The legacy `IDiracOperator` interface exposes private precomputed lookup tables (`getOmegaTable4()`) and breaks `const`-correctness by returning mutable matrix references from `const` member functions to allow MCMC mutation.

### 1.2 Goals
* **Research-Centric Workflow:** Provide a non-blocking **Stepper & Observer API** that allows researchers in C++ and Python to step through Markov sweeps, stream eigenvalue spectra, and measure integrated autocorrelation times $\tau_{\text{int}}$.
* **Modular Physics Engine:** Decouple the **Dirac State** (geometry), **Action** (energy & variations), and **Samplers** (MCMC / HMC algorithms).
* **Value Semantics & Memory Safety:** Transform `DiracOperator` into a clean, regular C++ value type (copyable, movable, strict `const` correctness, zero redundant heap wrappers).
* **Zero-Copy Python Interoperability:** Expose native Armadillo matrices and eigenvalue vectors directly to NumPy without memory copying.
* **Extensibility for Emerging Literature:** Ensure seamless addition of future physical models (e.g. Fermion Pfaffian effective actions, product geometries, HMC gradients).

### 1.3 Non-Goals
* Distributed MPI cluster scaling (single-node multi-threading via OpenMP is sufficient for current matrix sizes $N \le 256$).
* General-purpose symbolic algebra (RFL is focused purely on numerical spectral geometry).

---

## 2. Research Workflows & Scientific Requirements

### 2.1 Core Research Scenarios
1. **Scenario 1 (Interactive Exploration in Python):** A researcher sets up a spectral triple $(p,q,N)$ in a Jupyter notebook, runs 500 thermalization sweeps with automated dual-averaging step-size tuning, and plots the real-time eigenvalue density $\rho(\lambda)$.
2. **Scenario 2 (Automated Batch Sampling & Statistics):** An automated batch script executes 100,000 production sweeps across a parameter grid $(g_2, g_4)$, recording eigenvalue spectra every 10 sweeps to compute the spectral dimension $d_{\text{spec}}$ and edge eigenvalue statistics (e.g. comparing largest eigenvalue fluctuations against the [Tracy–Widom distribution](https://en.wikipedia.org/wiki/Tracy%E2%80%93Widom_distribution)).
3. **Scenario 3 (Extending Physics Actions):** A theorist implements a new action potential $S(D) = S_{\text{Barrett-Glaser}}(D) - \ln \text{Pf}(JD)$ by subclassing or providing a new Action policy without modifying any Monte Carlo sampler code.

### 2.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Physical & Mathematical Invariant |
| :--- | :--- | :--- |
| **REQ-001** | **Hermiticity & Signature Preservation** | Matrix variations $\delta M_k$ must strictly preserve $H_i^\dagger = H_i$ and $L_j^\dagger = -L_j$ dictated by signature signs $\epsilon_k \in \{+1, -1\}$. |
| **REQ-002** | **Exact Analytic Variation $\Delta S$** | Fast $\mathcal{O}(N^3)$ local variation must equal global recomputation $S(D + \delta D) - S(D)$ within floating-point tolerance ($< 10^{-10}$). |
| **REQ-003** | **Spectral Symmetry Preservation** | For symmetric geometries (where $\{\Gamma, D\} = 0$), the computed spectrum must satisfy $\{\lambda_i\} = \{-\lambda_i\}$ to machine precision ($< 10^{-12}$). |
| **REQ-004** | **Detailed Balance & Ergodicity** | The MCMC stepper must satisfy detailed balance and provide an integrated autocorrelation estimator $\tau_{\text{int}}$ to compute rigorous observable errors. |
| **REQ-005** | **Automated Step-Size Calibration** | Dual-averaging sweeps must tune the proposal scale to achieve target acceptance rates (e.g. $0.65 \pm 0.03$) during burn-in. |
| **REQ-006** | **Zero-Copy NumPy Interoperability** | C++ matrix and eigenvalue buffers must be exposed to Python/NumPy without memory copying or pointer slicing bugs. |

---

## 3. Architecture Decision Records (ADRs) & Trade-offs

### 3.1 ADR-1: Object Model & State Representation

| Criteria | Option A: Legacy Dynamic Interfaces (`IDiracOperator*`) | Option B: Regular Value Type (`DiracOperator`) | Option C: Opaque Handle / PIMPL |
| :--- | :--- | :--- | :--- |
| **Inlining in Hot Loops** | ❌ Poor (vtable indirection) | ✅ **Optimal (Direct inlining)** | ⚠️ Fair (pointer hop) |
| **Memory Layout** | ❌ Double heap (`unique_ptr<vector>`) | ✅ **Contiguous / Single allocation** | ⚠️ Extra heap allocation |
| **Const-Correctness** | ❌ Broken (`const` mutability hack) | ✅ **Strict C++ const-correctness** | ✅ Strict |
| **Ergonomics in C++** | ❌ Awkward (`std::move(unique_ptr)`) | ✅ **Natural (`DiracOperator D(p,q,N)`)** | ⚠️ Pointer semantics |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* In scientific computing, domain states are mathematical values. Treating `DiracOperator` as a regular value type allows natural copying, moving, stack allocation, and total compiler optimisation.

---

### 3.2 ADR-2: Simulation Execution & Orchestration

| Criteria | Option A: Monolithic Runner (`Simulation::run()`) | Option B: Stepper + Observer Pattern | Option C: Reactive / Coroutine Stream |
| :--- | :--- | :--- | :--- |
| **Caller Control** | ❌ Black box (all-or-nothing) | ✅ **Full control over loop & state** | ✅ Full control |
| **Observable Recording** | ❌ Hardcoded in simulation | ✅ **Decoupled via `ISimulationObserver`** | ⚠️ Complex setup |
| **Python Ergonomics** | ❌ Blocking GIL | ✅ **Pythonic generator (`iter_sweeps`)** | ⚠️ High binding complexity |
| **C++17 Compatibility** | ✅ Compatible | ✅ **Clean C++17 standard** | ❌ Requires C++20 coroutines |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* The Stepper + Observer pattern provides the ideal balance: researchers have direct control over their loops, while observers allow seamless, decoupled collection of eigenvalues, actions, and progress telemetry.

---

### 3.3 ADR-3: Action & Variation Architecture

| Criteria | Option A: Embedded in Sampler (Legacy) | Option B: Physics Action with Analytic Variations | Option C: Automatic Differentiation (AD) |
| :--- | :--- | :--- | :--- |
| **Separation of Concerns**| ❌ Physics tangled in MCMC | ✅ **Clean separation ($S, \Delta S, \nabla S$)** | ✅ Clean separation |
| **Computational Speed** | ✅ Fast analytical $\mathcal{O}(N^3)$ | ✅ **Fast analytical $\mathcal{O}(N^3)$** | ❌ 10–50x slower for matrix traces |
| **Reusability for HMC** | ❌ Impossible without duplication | ✅ **Directly reuses gradient $\nabla S$** | ✅ Automatic gradients |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* Analytical variation formulas are essential for MCMC performance. Encapsulating these formulas inside `BarrettGlaserAction` cleanly separates physical theory from numerical sampling and enables immediate reuse in HMC.

---

### 3.4 ADR-4: Action Extensibility for Fermionic Pfaffians & Multi-Term Models

In the literature (Barrett 2024), adding fermions modifies the partition function via the Grassmann integral:
$$Z = \int \mathcal{D}D \, e^{-S_B(D)} \text{Pfaffian}(JD) = \int \mathcal{D}D \, e^{-S_{\text{eff}}(D)}$$
where the acceptance probability factorizes as $\alpha = \min\left(1, e^{-\Delta S_B} \cdot \left| \frac{\text{Pfaffian}(JD')}{\text{Pfaffian}(JD)} \right|\right)$.

| Criteria | Option 1: Concrete Barrett-Glaser Core + Multiplicative Pfaffian Modifier | Option 2: Generic Composite Action Hierarchy (`IAction` + `add_term`) |
| :--- | :--- | :--- |
| **Execution Speed for Pure Barrett-Glaser** | ✅ Maximum (Direct inlining, zero virtual dispatch) | ⚠️ Lower (virtual dispatch overhead across terms) |
| **Numerical Stability for Fermions** | ✅ High (Computes Pfaffian ratio directly without log-divergences) | ⚠️ Moderate (Requires $-\ln|\text{Pf}|$ additive energy conversion) |
| **Mathematical Alignment** | ✅ Matches exact factorization of path integral $\Delta S_B + \text{ratio}$ | ⚠️ Treats non-local Pfaffians as standard polynomial traces |
| **Design Complexity** | ✅ Simple, direct, easily testable | ❌ High (Polymorphic action container, composite delta dispatch) |
| **Decision** | **Selected (Option 1)** | Deferred |

*Rationale:* Barrett-Glaser $S_B(D)$ is the permanent bosonic foundation of all current and future simulations. Keeping `BarrettGlaserAction` as a concrete, highly optimised class preserves maximum speed. When Fermionic simulations are introduced, a Pfaffian modifier is plugged in as a multiplicative acceptance weight, matching the exact mathematics of the Grassmann path integral.

---

## 4. Target Architecture & Component Design

*(Currently in collaborative discussion — component APIs and directory organization will be finalized next.)*
