# EP-1: Research-Driven Architecture Modernisation for RFL

* **Title:** Research-Driven Architecture Modernisation for RFL
* **Author:** Paul Druce
* **Status:** Accepted (In Progress)
* **Target Version:** RFL v0.3.0
* **Date:** 2026-08-30

---

## 1. Multi-Phase Implementation Tracker

This proposal spans the `v0.3.0` milestone.
The table below tracks the status of each implementation phase:

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Refactor `DiracOperator` into a regular value type with strict `const` correctness (ADR-1) | `v0.3.0` | [#12](https://github.com/pauldruce/RFL/issues/12) | ⏳ Scheduled |
| **Phase 2** | Extract analytic trace variations and gradients into dedicated `BarrettGlaserAction` (ADR-3) | `v0.3.0` | [#14](https://github.com/pauldruce/RFL/issues/14) | ⏳ Scheduled |
| **Phase 3** | Decouple simulation execution with non-blocking Stepper & Observer API and dual-averaging (ADR-2) | `v0.3.0` | [#13](https://github.com/pauldruce/RFL/issues/13) | ⏳ Scheduled |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Research Context
`RFL` was created to study **Finite Noncommutative Geometries (Random Fuzzy Spaces)** and **Spectral Triples** via Monte Carlo simulations. The core mathematical framework is governed by the Barrett-Glaser action (Barrett & Glaser 2016):

$$
S(D) = g_2 \mathrm{Tr}(D^2) + g_4 \mathrm{Tr}(D^4)
$$

where the Dirac operator is decomposed into Hermitian and anti-Hermitian matrices ($H_i, L_j$) coupled to Clifford gamma matrices:

$$
D = \sum_{i=1}^p \gamma^i \otimes H_i + \sum_{j=1}^q \gamma^{p+j} \otimes L_j
$$

The research programme expands into Fermion Functional Integrals (Barrett 2024), Product Geometries $S_F^2 \otimes \mathcal{F}$ (Barrett 2026), and Random Matrix Spectral Statistics.
Consequently, the legacy codebase exhibits structural bottlenecks:
1. **Opaque Execution Control:** The legacy `Simulation` class executes a monolithic, blocking loop. Researchers cannot inspect eigenvalue trajectories, stream observables to disk, or apply interactive stopping conditions in Jupyter notebooks.
2. **Entangled Physics & Sampling:** The analytical trace variation formulas $\Delta \mathrm{Tr}(D^2)$ and $\Delta \mathrm{Tr}(D^4)$ are embedded directly inside `Metropolis.cpp`. This prevents reusing the same physics for Hybrid Monte Carlo (HMC) or new action potentials. Researchers cannot easily add sextic terms $S_6 \propto \mathrm{Tr}(D^6)$ or fermion Pfaffians $-\ln \mathrm{Pf}(JD)$.
3. **Leaky Interface & Const Violation:** The legacy `IDiracOperator` interface exposes private precomputed lookup tables (`getOmegaTable4()`). It breaks `const`-correctness by returning mutable matrix references from `const` member functions to allow MCMC mutation.

### 2.2 Goals
* **Research-Centric Workflow:** Provide a non-blocking **Stepper & Observer API**. This interface allows researchers in C++ and Python to step through Markov sweeps, stream eigenvalue spectra, and measure autocorrelation times $\tau_{\text{int}}$.
* **Modular Physics Engine:** Decouple the **Dirac State** (geometry), **Action** (energy & variations), and **Samplers** (MCMC / HMC algorithms).
* **Value Semantics & Memory Safety:** Transform `DiracOperator` into a clean, regular C++ value type (copyable, movable, strict `const` correctness, zero redundant heap wrappers).
* **Zero-Copy Python Interoperability:** Expose native Armadillo matrices and eigenvalue vectors directly to NumPy without memory copying.
* **Extensibility for Emerging Literature:** Ensure seamless addition of future physical models (e.g. Fermion Pfaffian effective actions, product geometries, HMC gradients).

### 2.3 Non-Goals
* Distributed MPI cluster scaling (single-node multi-threading via OpenMP is sufficient for current matrix sizes $N \le 256$).
* General-purpose symbolic algebra (RFL is focused purely on numerical spectral geometry).

---

## 3. Research Workflows & Scientific Requirements

### 3.1 Core Research Scenarios
1. **Scenario 1 (Interactive Exploration in Python):** A researcher sets up a spectral triple $(p, q, N)$ in a Jupyter notebook. They run 500 thermalisation sweeps with automated dual-averaging step-size tuning, and plot the real-time eigenvalue density $\rho(\lambda)$.
2. **Scenario 2 (Automated Batch Sampling & Statistics):** An automated batch script executes 100,000 production sweeps across a parameter grid $(g_2, g_4)$. The script records eigenvalue spectra every 10 sweeps to compute the spectral dimension $d_{\text{spec}}$ and edge eigenvalue statistics (e.g. comparing largest eigenvalue fluctuations against the [Tracy–Widom distribution](https://en.wikipedia.org/wiki/Tracy%E2%80%93Widom_distribution)).
3. **Scenario 3 (Extending Physics Actions):** A theorist implements a new action potential $S(D) = S_{\text{BG}}(D) - \ln \mathrm{Pf}(JD)$. They achieve this by providing a new Action policy without modifying Monte Carlo sampler code.

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Lifecycle Status | Operational & Physical Invariant |
| :--- | :--- | :--- | :--- |
| **REQ-ARCH-001** | **Hermiticity & Signature Preservation** | Phase 1 (Implemented) | Matrix variations $\delta M_k$ must strictly preserve $H_i^\dagger = H_i$ and $L_j^\dagger = -L_j$ dictated by signature signs $\epsilon_k \in \lbrace +1, -1 \rbrace$. |
| **REQ-ARCH-002** | **Exact Analytic Variation** $\Delta S$ | Phase 2 (Planned) | Fast $\mathcal{O}(N^3)$ local variation must equal global recomputation $S(D + \delta D) - S(D)$ within scale-aware tolerance $\mathcal{O}(M \epsilon_{\mathrm{mach}} (\lvert S_f \rvert + \lvert S_i \rvert))$. |
| **REQ-ARCH-003** | **Spectral Symmetry Preservation** | Phase 1 (Implemented) | For symmetric geometries (where $\lbrace \Gamma, D \rbrace = 0$), the computed spectrum must satisfy $\lbrace \lambda_i \rbrace = \lbrace -\lambda_i \rbrace$ within scale-aware bound $\mathcal{O}(\epsilon_{\mathrm{mach}} \lVert D \rVert_2)$. |
| **REQ-ARCH-004** | **Detailed Balance & Ergodicity** | Phase 3 (Planned) | The MCMC stepper must satisfy detailed balance and provide an integrated autocorrelation estimator $\tau_{\text{int}}$ to compute rigorous observable errors. |
| **REQ-ARCH-005** | **Automated Step-Size Calibration** | Phase 3 (Planned) | Dual-averaging sweeps must tune the proposal scale to achieve target acceptance rates (e.g. $0.65 \pm 0.03$) during thermalisation. |
| **REQ-ARCH-006** | **Zero-Copy NumPy Interoperability** | Phase 1 (Implemented) | C++ matrix and eigenvalue buffers must be exposed to Python/NumPy without memory copying or pointer slicing bugs. |
| **REQ-ARCH-007** | **Active LTS Platform Policy** | Phase 1 (Implemented) | Maintain compatibility with active Ubuntu LTS (≥ 22.04), macOS (≥ 14), Armadillo (≥ 11.4), and conforming C++17 compilers. |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: Object Model & State Representation

| Criteria | Option A: Legacy Dynamic Interfaces (`IDiracOperator*`) | Option B: Regular Value Type (`DiracOperator`) | Option C: Opaque Handle / PIMPL |
| :--- | :--- | :--- | :--- |
| **Inlining in Hot Loops** | ❌ Poor (vtable indirection) | ✅ **Optimal (Direct inlining)** | ⚠️ Fair (pointer hop) |
| **Memory Layout** | ❌ Double heap (`unique_ptr<vector>`) | ✅ **Contiguous / Single allocation** | ⚠️ Extra heap allocation |
| **Const-Correctness** | ❌ Broken (`const` mutability hack) | ✅ **Strict C++ const-correctness** | ✅ Strict |
| **Ergonomics in C++** | ❌ Awkward (`std::move(unique_ptr)`) | ✅ **Natural (`DiracOperator D(p,q,N)`)** | ⚠️ Pointer semantics |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* In scientific computing, domain states are mathematical values. Treating `DiracOperator` as a regular value type allows natural copying, moving, stack allocation, and total compiler optimisation.

---

### 4.2 ADR-2: Simulation Execution & Orchestration

| Criteria | Option A: Monolithic Runner (`Simulation::run()`) | Option B: Stepper + Observer Pattern | Option C: Reactive / Coroutine Stream |
| :--- | :--- | :--- | :--- |
| **Caller Control** | ❌ Black box (all-or-nothing) | ✅ **Full control over loop & state** | ✅ Full control |
| **Observable Recording** | ❌ Hardcoded in simulation | ✅ **Decoupled via `ISimulationObserver`** | ⚠️ Complex setup |
| **Python Ergonomics** | ❌ Blocking GIL | ✅ **Pythonic generator (`iter_sweeps`)** | ⚠️ High binding complexity |
| **C++17 Compatibility** | ✅ Compatible | ✅ **Clean C++17 standard** | ❌ Requires C++20 coroutines |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* The Stepper + Observer pattern provides the ideal balance. Researchers have direct control over their loops, while observers allow decoupled collection of eigenvalues, actions, and progress telemetry.

---

### 4.3 ADR-3: Action & Variation Architecture

| Criteria | Option A: Embedded in Sampler (Legacy) | Option B: Physics Action with Analytic Variations | Option C: Automatic Differentiation (AD) |
| :--- | :--- | :--- | :--- |
| **Separation of Concerns**| ❌ Physics tangled in MCMC | ✅ **Clean separation** ($S, \Delta S, \nabla S$) | ✅ Clean separation |
| **Computational Speed** | ✅ Fast analytical $\mathcal{O}(N^3)$ | ✅ **Fast analytical** $\mathcal{O}(N^3)$ | ❌ 10–50x slower for matrix traces |
| **Reusability for HMC** | ❌ Impossible without duplication | ✅ **Directly reuses gradient** $\nabla S$ | ✅ Automatic gradients |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* Analytical variation formulas are essential for MCMC performance. Encapsulating these formulas inside `BarrettGlaserAction` cleanly separates physical theory from numerical sampling and enables immediate reuse in HMC.

---

### 4.4 ADR-4: Action Extensibility for Fermionic Pfaffians & Multi-Term Models

In the literature (Barrett 2024), adding fermions modifies the partition function via the Grassmann integral:

$$
Z = \int \mathcal{D}D \thinspace \mathrm{e}^{-S_B(D)} \mathrm{Pfaffian}(JD) = \int \mathcal{D}D \thinspace \mathrm{e}^{-S_{\mathrm{eff}}(D)}
$$

where the acceptance probability factorises as:

$$
\alpha = \min\left(1, \mathrm{e}^{-\Delta S_B} \cdot \left| \frac{\mathrm{Pfaffian}(JD')}{\mathrm{Pfaffian}(JD)} \right|\right)
$$

| Criteria | Option 1: Concrete Barrett-Glaser Core + Multiplicative Pfaffian Modifier | Option 2: Generic Composite Action Hierarchy (`IAction` + `add_term`) |
| :--- | :--- | :--- |
| **Execution Speed for Pure Barrett-Glaser** | ✅ Maximum (Direct inlining, zero virtual dispatch) | ⚠️ Lower (virtual dispatch overhead across terms) |
| **Numerical Stability for Fermions** | ✅ High (Computes Pfaffian ratio directly without log-divergences) | ⚠️ Moderate (Requires $-\ln\lvert\text{Pf}\rvert$ additive energy conversion) |
| **Mathematical Alignment** | ✅ Matches exact factorisation of path integral $\Delta S_B + \text{ratio}$ | ⚠️ Treats non-local Pfaffians as standard polynomial traces |
| **Design Complexity** | ✅ Simple, direct, easily testable | ❌ High (Polymorphic action container, composite delta dispatch) |
| **Decision** | **Selected (Option 1)** | Deferred |

*Rationale:* Barrett-Glaser $S_B(D)$ is the permanent bosonic foundation of all current and future simulations. Keeping `BarrettGlaserAction` as a concrete, highly optimised class preserves maximum speed. When Fermionic simulations are introduced, a Pfaffian modifier is plugged in as a multiplicative acceptance weight, matching the exact mathematics of the Grassmann path integral.

---

## 5. Target Architecture & Component Design

*(Currently in collaborative discussion — component APIs and directory organisation will be finalised next.)*
