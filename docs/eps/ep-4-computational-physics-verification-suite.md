# EP-4: Computational Physics Verification Suite & Scientific Validation Standards

<!-- cspell:words Bonferroni GROMACS PRNG Smirnov Zenodo conftest unitaries -->

* **Title:** Computational Physics Verification Suite & Scientific Validation Standards
* **Author:** Paul Druce
* **Status:** Draft
* **Target Versions:** RFL v0.3.0 (Phase 1), v0.4.0 (Phase 2)
* **Date:** 2026-09-17

---

## 1. Multi-Phase Implementation Tracker

This proposal establishes a scientific verification test suite for RFL.
The table below tracks the status of each implementation phase:

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Port `delta24` action tests; implement complete 8-fold KO spectral axioms; add unitary gauge invariance; add Gaussian limit checks. | `v0.3.0` | ⏳ Scheduled | 💡 Draft |
| **Phase 2** | Add 1-matrix Riemann-Hilbert tests; detailed balance flux checks; Simulation-Based Calibration; rank-normalized split-$\hat{R}$; Zenodo DOI. | `v0.4.0` | 💡 Planned | 💡 Draft |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Research Context
Computational physics software requires rigorous verification against exact mathematical theorems and physical limits.
Standard unit tests only verify that code runs without runtime errors or crashes.
However, unit tests do not prove mathematical or physical correctness.
Markov Chain Monte Carlo (MCMC) algorithms can produce plausible numbers while silently violating geometric axioms or Boltzmann statistics.

Random Noncommutative Geometry studies quantum spaces through path integrals over Dirac operators.
The partition function integrates over the space of Dirac operators in a finite spectral triple:
$$Z = \int \mathcal{D}D \, \mathrm{e}^{-S(D)}$$
The Barrett-Glaser spectral action governs this ensemble:
$$S(D) = g_2 \mathrm{Tr}(D^2) + g_4 \mathrm{Tr}(D^4)$$

To establish scientific validity, RFL must verify six fundamental pillars:
1. **Axiomatic Geometric Invariants:** The assembled Dirac operator must satisfy all spectral triple axioms across all 8 KO dimensions. These axioms include Hermiticity, real structure relations ($J^2 = \epsilon', DJ = \epsilon JD, J\Gamma = \epsilon'' \Gamma J$), and chirality grading ($\{\Gamma, D\} = 0$).
2. **Continuous Symmetries & Gauge Invariance:** The spectral action and eigenvalue spectra must be invariant under unitary transformations $D \to UDU^\dagger$ for $U \in \mathrm{U}(N)$.
3. **Internal Energy & Derivative Invariants:** Incremental action updates $\Delta S$ (`delta24`) must match full action evaluations $S(D_f) - S(D_i)$ to double precision. Analytic variations must match central finite differences to $\mathcal{O}(h^2)$.
4. **Exact Limiting Theorems:** In solvable limits, simulated observables must reproduce analytical predictions. In the Gaussian regime ($g_4 = 0$), spectral moments must match the self-convolution of the Wigner semicircle distribution. One-matrix reductions ($(0, 1)$) must match exact Riemann-Hilbert solutions.
5. **Statistical Mechanics & Ergodicity:** Markov chain updates must satisfy detailed balance with respect to the Boltzmann weight. Multi-chain simulations must verify convergence using rank-normalized folded split-$\hat{R}$ and Effective Sample Size diagnostics ($\mathrm{ESS} \ge 400$).
6. **Metamorphic Invariants for Unsolved Regimes:** For general interactive potentials where exact analytical solutions do not exist, simulations must satisfy parameter rescaling and coupling monotonicity relations.

Automating these checks ensures that RFL provides reliable, mathematically sound foundations for research.
Researchers and peer reviewers can independently confirm that the software samples the true physical distribution.

### 2.2 Goals
* **Goal 1:** Verify the incremental action update $\Delta S$ (`delta24`) against full action evaluation to double precision ($< 10^{-10}$).
* **Goal 2:** Verify fundamental spectral triple axioms (Hermiticity, complete 8-fold KO real structure table, chirality) for all generated configurations.
* **Goal 3:** Verify unitary gauge invariance ($U(N)$ symmetry) for action evaluations and eigenvalue spectra under random Haar unitary transformations.
* **Goal 4:** Verify analytical action derivatives and variations against numerical central finite differences to $\mathcal{O}(h^2)$.
* **Goal 5:** Implement an automated Gaussian-limit test comparing MCMC eigenvalue moments to the exact Wigner semicircle self-convolution and Simulation-Based Calibration (SBC).
* **Goal 6:** Validate 1-matrix models ($(0, 1)$) against exact analytical Riemann-Hilbert solutions from published literature.
* **Goal 7:** Provide automated MCMC statistical diagnostics including integrated autocorrelation time $\tau_{\mathrm{int}}$ and rank-normalized folded split-$\hat{R}$ with $\mathrm{ESS}_{\mathrm{bulk}} \ge 400$.

### 2.3 Non-Goals
* Re-implementing external statistical packages inside RFL C++ core.
* Modifying existing Clifford algebra matrix representations.
* Adding new interaction potentials beyond quartic polynomials in this proposal.

---

## 3. Research Workflows & Scientific Requirements

### 3.1 Core Research Scenarios
1. **Scenario 1 (Continuous Integration Verification):**
   A developer modifies C++ matrix algebra.
   The fast test runner executes deterministic algebraic and geometric invariant tests.
   The runner confirms that action calculations, unitary gauge symmetry, and Clifford relations match exact mathematical identities within machine precision.

2. **Scenario 2 (Academic Peer-Review Audit):**
   A journal reviewer downloads RFL from GitHub.
   The reviewer executes `pytest tests/physics/` or `ctest`.
   The suite outputs statistical verification metrics and confirms detailed balance without requiring manual configuration.

3. **Scenario 3 (Metamorphic Testing of Interactive Ensembles):**
   A researcher executes simulations in parameter regimes with no analytical solution ($g_4 > 0, N > 1$).
   The test suite checks parameter rescaling invariance and coupling monotonicity to ensure physical behaviour.

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Physical & Mathematical Invariant | Target Tier |
| :--- | :--- | :--- | :--- |
| **REQ-001** | **Incremental Action Invariant** | $|(S_f - S_i) - \Delta S| < 10^{-10}$ for all single-element matrix updates. | Tier 1 (Smoke) |
| **REQ-002** | **Dirac Hermiticity** | $\|D - D^\dagger\|_F < 10^{-14}$ for assembled Dirac operators. | Tier 1 (Smoke) |
| **REQ-003** | **Complete KO Real Structure** | $\|J^2 - \epsilon' I\|_F < 10^{-14}$, $\|D J - \epsilon J D\|_F < 10^{-14}$, and $\|J \Gamma - \epsilon'' \Gamma J\|_F < 10^{-14}$ across all 8 KO dimensions. | Tier 1 (Smoke) |
| **REQ-004** | **Chirality & Spectral Anti-Symmetry** | $\|\Gamma D + D \Gamma\|_F < 10^{-14}$ and $\{\lambda_i\} \equiv \{-\lambda_i\}$ to machine precision for even spectral triples. | Tier 1 (Smoke) |
| **REQ-005** | **Gaussian Limit Moments** | MCMC eigenvalue moments match Wigner self-convolution within $3$ standard errors. | Tier 2 (Integration) |
| **REQ-006** | **Detailed Balance** | Transition probabilities satisfy microscopic reversibility and coarse-grained state-flux balance $N_{A \to B} \approx N_{B \to A}$. | Tier 2 (Integration) |
| **REQ-007** | **Chain Convergence** | Rank-normalized folded split-$\hat{R} < 1.05$ with $\mathrm{ESS}_{\mathrm{bulk}} \ge 400$ across independent chains. | Tier 3 (Validation) |
| **REQ-008** | **Unitary Gauge Invariance** | $\|S(U D U^\dagger) - S(D)\| < 10^{-12}$ and $\max_i \|\lambda_i(U D U^\dagger) - \lambda_i(D)\| < 10^{-12}$ for Haar unitary $U \in \mathrm{U}(N)$. | Tier 1 (Smoke) |
| **REQ-009** | **Action Derivatives** | Analytical matrix variations match central finite differences to $\mathcal{O}(h^2)$ across all matrix elements. | Tier 1 (Smoke) |
| **REQ-010** | **Simulation-Based Calibration** | Posterior rank statistics from Gaussian ensemble MCMC pass Kolmogorov-Smirnov uniformity test ($p > 0.01$). | Tier 3 (Validation) |
| **REQ-011** | **Metamorphic Invariants** | Rescaling $(g_2, g_4) \to (\alpha^{-2} g_2, \alpha^{-4} g_4)$ and monotonicity of $\langle \mathrm{Tr}(D^2) \rangle$ hold for general potentials. | Tier 2 (Integration) |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: Test Suite Architecture (C++ vs Python)

| Criteria | Option A: Pure C++ GoogleTest | Option B: Two-Tier (C++ Core + Python Pytest) (Selected) |
| :--- | :--- | :--- |
| **Statistical Analysis** | High development cost in C++ | Native access to SciPy and NumPy |
| **Execution Speed** | Fast for micro-benchmarks | Fast C++ execution via `pyrfl` bindings |
| **Analytical Comparison** | Difficult to express integrals | Simple quadrature and curve fitting |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:*
C++ tests in `src/core/tests/` verify deterministic algebraic identities and geometric invariants.
Python tests in `tests/physics/` verify statistical distributions, moments, and Riemann-Hilbert curves using SciPy.

### 4.2 ADR-2: Three-Tier Testing Pyramid & Stochastic Flakiness Policy

| Criteria | Option A: Monolithic CI Test Suite | Option B: Three-Tier Testing Pyramid (Selected) |
| :--- | :--- | :--- |
| **CI Feedback Speed** | Slow (> 30 mins on every PR) | Fast (< 2 mins for PR smoke gate) |
| **Flakiness Management** | High (frequent stochastic false alarms) | Controlled via Holm-Bonferroni correction and two-stage re-seed protocol |
| **Statistical Rigour** | Weak (compromised to keep CI fast) | High (extended sampling on nightly runners) |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:*
Following established practices in GROMACS and Stan:
* **Tier 1 (Smoke / Deterministic Gate, < 2 minutes):** Executes on every commit and PR. Uses deterministic PRNG seeds to test algebraic identities, Hermiticity, gauge invariance, and finite differences.
* **Tier 2 (Physics Integration Gate, ~3–5 minutes):** Executes on PR merge to `main`. Uses fixed seeds to test short MCMC chains, loose-tolerance Gaussian moments ($4\sigma$), and state-flux balance.
* **Tier 3 (Scientific Validation Suite, Nightly / Release Gate, ~30–60 minutes):** Executes high-statistics runs, multi-chain split-$\hat{R}$, Riemann-Hilbert curve fits, and Simulation-Based Calibration.
* **Flakiness Control Protocol:** Stochastic tests apply the Holm-Bonferroni correction to prevent false discovery inflation. If a stochastic test fails ($p < \alpha$), CI triggers an automated second run with an independent seed and double chain length before failing.

---

## 5. Target Architecture & Component Contracts

### 5.1 C++ Core Verification (`src/core/tests/`)
C++ tests focus on deterministic invariants and algebraic correctness:
* `tDelta.cpp`: Verifies incremental action update `delta24` against full action difference across all Clifford types $(p, q)$ with $p+q \le 4$.
* `tAxioms.cpp`: Verifies Hermiticity ($D = D^\dagger$), complete 8-fold KO table ($J^2, JD, J\Gamma$), and chirality ($\{\Gamma, D\} = 0$).
* `tGaugeInvariance.cpp`: Verifies that random Haar unitary transformations $D \to UDU^\dagger$ preserve action values and eigenvalue spectra.
* `tDerivatives.cpp`: Verifies analytical matrix variations against central finite differences.

### 5.2 Python Physics Suite (`tests/physics/`)
Python tests focus on statistical physics, limiting theorems, and asymptotic distributions:
```
tests/physics/
├── __init__.py
├── conftest.py
├── test_axiomatic_invariants.py
├── test_gauge_invariance.py
├── test_derivative_consistency.py
├── test_gaussian_limit.py
├── test_gaussian_sbc.py
├── test_riemann_hilbert_1matrix.py
├── test_detailed_balance.py
└── test_metamorphic_invariants.py
```

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Test Description | Target Invariant | Target Tier |
| :--- | :--- | :--- | :--- |
| **Algebraic Invariants** | `ctest --test-dir build -R DeltaTests` | Verifies $\Delta S$ matches full action differences. | Tier 1 (PR Gate) |
| **Spectral Axioms** | `pytest tests/physics/test_axiomatic_invariants.py` | Verifies Hermiticity, KO table, and chirality. | Tier 1 (PR Gate) |
| **Gauge Invariance** | `pytest tests/physics/test_gauge_invariance.py` | Verifies $S(UDU^\dagger) = S(D)$ and spectral invariance. | Tier 1 (PR Gate) |
| **Derivative Consistency** | `pytest tests/physics/test_derivative_consistency.py` | Verifies variations match central finite differences. | Tier 1 (PR Gate) |
| **Gaussian Limit** | `pytest tests/physics/test_gaussian_limit.py` | Verifies convergence to Wigner convolution law. | Tier 2 (Merge Gate) |
| **Detailed Balance** | `pytest tests/physics/test_detailed_balance.py` | Verifies microscopic reversibility and state-flux balance. | Tier 2 (Merge Gate) |
| **Metamorphic Relations** | `pytest tests/physics/test_metamorphic_invariants.py` | Verifies parameter rescaling and coupling monotonicity. | Tier 2 (Merge Gate) |
| **1-Matrix Benchmark** | `pytest tests/physics/test_riemann_hilbert_1matrix.py` | Verifies agreement with Riemann-Hilbert analytical solutions. | Tier 3 (Nightly) |
| **MCMC Calibration (SBC)** | `pytest tests/physics/test_gaussian_sbc.py` | Verifies uniform posterior ranks in Gaussian ensembles. | Tier 3 (Nightly) |
| **Chain Convergence** | `pytest tests/physics/test_chain_convergence.py` | Verifies rank-normalized split-$\hat{R} < 1.05$ and $\mathrm{ESS} \ge 400$. | Tier 3 (Nightly) |

---

## 7. Phased Delivery Plan

### Phase 1: Core Algebraic & Limiting Physics Verification
* **Target Version:** `v0.3.0`
* **Alignment:** Implemented alongside [EP-1](ep-1-core-architecture-modernisation.md) architecture modernisation.
* **Tasks:**
  1. Create `src/core/tests/tDelta.cpp` testing `delta24` across all Clifford types $(p, q)$ with $p+q \le 4$.
  2. Implement `tests/physics/test_axiomatic_invariants.py` asserting $D = D^\dagger$, the complete 8-fold KO table, and $\{\Gamma, D\} = 0$.
  3. Implement `tests/physics/test_gauge_invariance.py` asserting unitary invariance $S(UDU^\dagger) = S(D)$ under random Haar unitaries.
  4. Implement `tests/physics/test_derivative_consistency.py` testing variations against central finite differences.
  5. Implement `tests/physics/test_gaussian_limit.py` asserting Gaussian trace moments match Wigner self-convolution within statistical tolerance.
  6. Add Tier 1 and Tier 2 test suites to `.github/workflows/ci.yml`.

### Phase 2: Statistical Rigour, Benchmarks & Open Science
* **Target Version:** `v0.4.0`
* **Tasks:**
  1. Implement `tests/physics/test_riemann_hilbert_1matrix.py` verifying $(0, 1)$ spectral densities against exact analytical curves.
  2. Implement `tests/physics/test_detailed_balance.py` testing microscopic reversibility and coarse-grained state-flux balance.
  3. Implement `tests/physics/test_gaussian_sbc.py` validating the MCMC sampler with Simulation-Based Calibration.
  4. Implement rank-normalized folded split-$\hat{R}$ and $\mathrm{ESS}_{\mathrm{bulk}}$ diagnostics in Python analysis utilities.
  5. Implement `tests/physics/test_metamorphic_invariants.py` testing parameter rescaling and coupling monotonicity.
  6. Create benchmark scripts to recreate published matrix model scaling curves and archive golden reference datasets.
  7. Configure automated Zenodo DOI archiving upon GitHub release tags.
