# EP-4: Computational Physics Verification Suite & Scientific Validation Standards

<!-- cspell:words Bonferroni CmdStan Feynman GROMACS Hartree Hellmann Higham LAMMPS Liouville MILC OpenMM PRNG Psi4 PySCF Smirnov USQCD Wielandt Wilkinson Zenodo autodiff conftest microcanonical plaquette plaquettes pseudofermion significand symplectic unitaries -->

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
| **Phase 1** | Port `delta24` action tests; implement complete 8-fold KO spectral axioms; add unitary gauge invariance; add Gaussian limit checks. | `v0.3.0` | [#64](https://github.com/pauldruce/RFL/pull/64) | 💡 Draft |
| **Phase 2** | Add 1-matrix Riemann-Hilbert tests; detailed balance flux checks; Simulation-Based Calibration; rank-normalized split $\hat{R}$; Zenodo DOI. | `v0.4.0` | 💡 Planned | 💡 Draft |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Research Context
Computational physics software requires verification across both deterministic mathematics and statistical physics.
Existing unit tests in RFL successfully verify local algebraic relations, Clifford anticommutator identities, and matrix symmetries.
However, isolated unit tests cannot verify collective, asymptotic, or stochastic physical behaviour.
Markov Chain Monte Carlo (MCMC) samplers can satisfy local algebraic checks while silently violating Boltzmann statistics, ergodicity, or limiting spectral theorems.

Random Noncommutative Geometry studies quantum spaces through path integrals over Dirac operators.
The partition function integrates over the space of Dirac operators in a finite spectral triple:


$$
Z = \int \mathcal{D}D \, \mathrm{e}^{-S(D)}
$$


The Barrett-Glaser spectral action governs this ensemble:


$$
S(D) = g_2 \mathrm{Tr}(D^2) + g_4 \mathrm{Tr}(D^4)
$$



To establish scientific validity, RFL must verify six fundamental pillars:
1. **Axiomatic Geometric Invariants:** The assembled Dirac operator must satisfy all spectral triple axioms across all 8 KO dimensions. These axioms include Hermiticity, real structure relations ($J^2 = \epsilon', DJ = \epsilon JD, J\Gamma = \epsilon'' \Gamma J$), and chirality grading ($\{\Gamma, D\} = 0$).
2. **Continuous Symmetries & Gauge Invariance:** The spectral action and eigenvalue spectra must be invariant under unitary transformations $D \to UDU^\dagger$ for $U \in \mathrm{U}(N)$.
3. **Internal Energy & Derivative Invariants:** Incremental action updates $\Delta S$ (`delta24`) must match full action evaluations $S(D_f) - S(D_i)$ within scale-aware cancellation bounds. Analytic variations must match central finite differences to optimal floating-point precision $\mathcal{O}(\epsilon_{\mathrm{mach}}^{2/3}) \lVert \nabla S \rVert$. In MCMC and HMC, these invariants guarantee that the sampler experiences true physical potential energy and conservative forces rather than unphysical numerical drift.
4. **Exact Limiting Theorems:** In solvable limits, simulated observables must reproduce analytical predictions. In the Gaussian regime with $g_4 = 0$, spectral moments must match the self-convolution of the Wigner semicircle distribution. One-matrix reductions of signature $(0, 1)$ must match exact Riemann-Hilbert solutions.
5. **Statistical Mechanics & Ergodicity:** Markov chain updates must satisfy detailed balance with respect to the Boltzmann weight. Multi-chain simulations must verify convergence using rank-normalized folded split $\hat{R}$ and Effective Sample Size diagnostics ($\mathrm{ESS} \ge 400$).
6. **Physical Scaling Invariants for General Regimes:** In interactive regimes without analytical solutions, simulations must satisfy exact parameter rescaling and coupling monotonicity laws.

#### The Unit-Test-First Principle for Physical Invariants
Fast unit tests provide immediate deterministic feedback with zero stochastic flakiness.
Whenever a physical invariant can be verified on a single configuration or state transition, developers must write a unit test.
Single-step action updates, geometric axioms, gauge transformations, and finite-difference derivatives must execute as unit tests.
Longer stochastic simulations must only test collective properties, ergodicity, and asymptotic distributions that unit tests cannot cover.

Automating these checks ensures that RFL provides reliable, mathematically sound foundations for research.
Researchers and peer reviewers can independently confirm that the software samples the true physical distribution.

### 2.2 Goals
* **Goal 1:** Verify the incremental action update $\Delta S$ (`delta24`) against full action evaluation to scale-aware precision.
* **Goal 2:** Verify fundamental spectral triple axioms (Hermiticity, complete 8-fold KO real structure table, chirality) for all generated configurations.
* **Goal 3:** Verify unitary gauge invariance ($U(N)$ symmetry) for action evaluations and eigenvalue spectra under random Haar unitary transformations.
* **Goal 4:** Verify analytical action derivatives and variations against numerical central finite differences to $\mathcal{O}(\epsilon_{\mathrm{mach}}^{2/3}) \lVert \nabla S \rVert$.
* **Goal 5:** Implement an automated Gaussian-limit test comparing MCMC eigenvalue moments to the exact Wigner semicircle self-convolution and Simulation-Based Calibration (SBC).
* **Goal 6:** Validate 1-matrix models of signature $(0, 1)$ against exact analytical Riemann-Hilbert solutions from published literature.
* **Goal 7:** Provide automated MCMC statistical diagnostics including integrated autocorrelation time $\tau_{\mathrm{int}}$ and rank-normalized folded split $\hat{R}$ with $\mathrm{ESS}_{\mathrm{bulk}} \ge 400$.

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

3. **Scenario 3 (Verification in General Interactive Regimes):**
   A researcher runs simulations in parameter regimes without analytical solutions with $g_4 > 0$ and $N > 1$.
   The test suite verifies parameter rescaling invariance and coupling monotonicity to confirm physical behaviour.

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Physical & Mathematical Invariant | Target Tier |
| :--- | :--- | :--- | :--- |
| **REQ-001** | **Incremental Action Invariant** | $\lvert (S_f - S_i) - \Delta S \rvert \le c_1 M \epsilon_{\mathrm{mach}} (\lvert S_f \rvert + \lvert S_i \rvert) + 10^{-12}$ (Section 3.3). | Tier 1 (Smoke) |
| **REQ-002** | **Dirac Hermiticity** | $\lVert D - D^\dagger \rVert_F \le c_2 (p+q) \epsilon_{\mathrm{mach}} \lVert D \rVert_F$ for assembled Dirac operators. | Tier 1 (Smoke) |
| **REQ-003** | **Complete KO Real Structure** | $\lVert J^2 - \epsilon' I \rVert_F \le c_3 \epsilon_{\mathrm{mach}}$, $\lVert D J - \epsilon J D \rVert_F \le c_3 \epsilon_{\mathrm{mach}} \lVert D \rVert_F$, and $\lVert J \Gamma - \epsilon'' \Gamma J \rVert_F \le c_3 \epsilon_{\mathrm{mach}}$ across all 8 KO dimensions. | Tier 1 (Smoke) |
| **REQ-004** | **Chirality & Spectral Anti-Symmetry** | $\lVert \Gamma D + D \Gamma \rVert_F \le 2 \epsilon_{\mathrm{mach}} \lVert D \rVert_F$ and $\{\lambda_i\} \equiv \{-\lambda_i\}$ to machine precision for even spectral triples. | Tier 1 (Smoke) |
| **REQ-005** | **Gaussian Limit Moments** | MCMC eigenvalue moments match Wigner self-convolution within $3$ standard errors. | Tier 2 (Integration) |
| **REQ-006** | **Detailed Balance** | Transition probabilities satisfy microscopic reversibility and coarse-grained state-flux balance $N_{A \to B} \approx N_{B \to A}$. | Tier 2 (Integration) |
| **REQ-007** | **Chain Convergence** | Rank-normalized folded split $\hat{R} < 1.05$ with $\mathrm{ESS}_{\mathrm{bulk}} \ge 400$ across independent chains. | Tier 3 (Validation) |
| **REQ-008** | **Unitary Gauge Invariance** | $\lvert S(U D U^\dagger) - S(D) \rvert \le c_4 M^2 \epsilon_{\mathrm{mach}} \lvert S(D) \rvert$ and $\max_i \lvert \lambda_i(U D U^\dagger) - \lambda_i(D) \rvert \le c_5 M \epsilon_{\mathrm{mach}} \lVert D \rVert_2$ for Haar unitary $U \in \mathrm{U}(N)$. | Tier 1 (Smoke) |
| **REQ-009** | **Action Derivatives** | Analytical matrix variations match central finite differences to $\mathcal{O}(\epsilon_{\mathrm{mach}}^{2/3}) \lVert \nabla S \rVert$ (Section 3.3). | Tier 1 (Smoke) |
| **REQ-010** | **Simulation-Based Calibration** | Posterior rank statistics from Gaussian ensemble MCMC pass Kolmogorov-Smirnov uniformity test with $p > 0.01$. | Tier 3 (Validation) |
| **REQ-011** | **Physical Scaling Laws** | Rescaling $(g_2, g_4) \to (\alpha^{-2} g_2, \alpha^{-4} g_4)$ and monotonicity of $\langle \mathrm{Tr}(D^2) \rangle$ hold for general potentials. | Tier 2 (Integration) |

---

### 3.3 Numerical Error Bounds & Precision Justifications

Hardcoded absolute tolerances fail when matrix scale or dimensions change.
Following LAPACK and Higham (2002), RFL grounds all tolerances in IEEE 754 double precision machine epsilon:


$$
\epsilon_{\mathrm{mach}} = 2^{-52} \approx 2.2204 \times 10^{-16}
$$


and the Hilbert space dimension of the Dirac operator:


$$
M = \dim(D) = 2^{\lfloor (p+q)/2 \rfloor} N
$$



#### 1. Structural Algebraic Invariants (REQ-002, REQ-003, REQ-004)
The Dirac operator is assembled by summing $p+q$ Kronecker tensor products:


$$
D = \sum_{k=1}^{p+q} \gamma^k \otimes M_k
$$


Because $\gamma^k$ contains exact matrix elements in $\{0, \pm 1, \pm i\}$, errors arise only from floating-point additions.
Summing $p+q \le 8$ terms accumulates at most $\mathcal{O}((p+q)\epsilon_{\mathrm{mach}})$ rounding error:


$$
\frac{\lVert D - D^\dagger \rVert_F}{\lVert D \rVert_F} \le c_2 (p+q) \epsilon_{\mathrm{mach}} \approx 10^{-14}
$$


Similarly, charge conjugation $J$ and grading $\Gamma$ represent signed index permutations.
Matrix multiplication by these operators incurs zero cancellation and minimal floating-point error bounded by $c_3 \epsilon_{\mathrm{mach}} \|D\|_F$.

#### 2. Catastrophic Cancellation in Action Differences (REQ-001)
When updating a single matrix element, the full action evaluation computes:


$$
\Delta S_{\mathrm{full}} = S(D_f) - S(D_i)
$$


The action values $S(D_f)$ and $S(D_i)$ can exceed $10^4$.
However, the single-element variation $\Delta S$ can be small ($\approx 10^{-2}$).
Subtracting two large, nearly equal floating-point numbers causes catastrophic cancellation.
The calculation loses $\log_{10}(|S| / |\Delta S|) \approx 6$ decimal digits of precision.
In contrast, `delta24` computes $\Delta S$ directly from local matrix variations, avoiding subtraction.
The discrepancy between `delta24` and $\Delta S_{\mathrm{full}}$ is bounded by the cancellation error of the full traces:


$$
\lvert (S_f - S_i) - \Delta S \rvert \le c_1 M \epsilon_{\mathrm{mach}} (\lvert S_f \rvert + \lvert S_i \rvert) + 10^{-12}
$$


For $M \approx 100$ and $|S| \approx 10^4$, this error is approximately $10^{-10}$ to $10^{-9}$.

#### 3. Unitary Gauge & Spectral Invariance (REQ-008)
Unitary transformation $D \to U D U^\dagger$ involves matrix multiplication.
Numerical matrix multiplication of dimension $M$ incurs rounding error bounded by Wilkinson's theorem:


$$
\lVert fl(UDU^\dagger) - UDU^\dagger \rVert_F \le c M \epsilon_{\mathrm{mach}} \lVert D \rVert_F
$$


Trace evaluation of $(UDU^\dagger)^4$ accumulates additional error of order $M^2 \epsilon_{\mathrm{mach}} \lvert S(D) \rvert$.
For eigenvalue spectra, the Hoffman-Wielandt and Weyl perturbation theorems bound eigenvalue drift by the operator 2-norm:


$$
\lvert \lambda_i(U D U^\dagger) - \lambda_i(D) \rvert \le \lVert fl(UDU^\dagger) - UDU^\dagger \rVert_2 \le c_5 M \epsilon_{\mathrm{mach}} \lVert D \rVert_2
$$


For $M \le 100$ and $\|D\|_2 \approx 10$, eigenvalue drift remains below $10^{-13}$.

#### 4. Theoretical Bound for Central Finite Differences (REQ-009)
Central finite difference approximations balance truncation error and floating-point cancellation:


$$
\frac{\partial S}{\partial M_{IJ}} \approx \frac{S(M_{IJ} + h) - S(M_{IJ} - h)}{2h}
$$


The total error is bounded by:


$$
\mathrm{Error}(h) \le \frac{h^2}{6} \lvert S''' \rvert + \frac{\epsilon_{\mathrm{mach}}}{h} \lvert S \rvert
$$


Minimising this error yields the optimal step size:


$$
h^* = \left( \frac{3 \epsilon_{\mathrm{mach}} \lvert S \rvert}{\lvert S''' \rvert} \right)^{1/3} \approx \mathcal{O}(\epsilon_{\mathrm{mach}}^{1/3}) \approx 6 \times 10^{-6}
$$


Evaluating at $h^*$ yields the minimum achievable error:


$$
\mathrm{Error}_{\mathrm{min}} \approx \mathcal{O}(\epsilon_{\mathrm{mach}}^{2/3}) \lvert S \rvert \approx (2.22 \times 10^{-16})^{2/3} \lvert S \rvert \approx 3.6 \times 10^{-11} \lvert S \rvert
$$


Therefore, asserting tolerances tighter than $10^{-10}$ for finite differences in double precision is mathematically invalid.
Tests must use a relative tolerance scaled by $\mathcal{O}(\epsilon_{\mathrm{mach}}^{2/3})$.

#### 5. Summary of Numerical Bounds

| Requirement ID | Error Bound Mechanism | Mathematical Error Bound | Empirical Constant |
| :--- | :--- | :--- | :--- |
| **REQ-001** | Trace Cancellation | $\le c_1 M \epsilon_{\mathrm{mach}} (\lvert S_f \rvert + \lvert S_i \rvert) + 10^{-12}$ | $c_1 \approx 10$ |
| **REQ-002** | Additive Hermiticity | $\le c_2 (p+q) \epsilon_{\mathrm{mach}} \lVert D \rVert_F$ | $c_2 \approx 5$ |
| **REQ-003** | Permutation Invariance | $\le c_3 \epsilon_{\mathrm{mach}} \lVert D \rVert_F$ | $c_3 \approx 5$ |
| **REQ-004** | Sign Anti-commutation | $\le 2 \epsilon_{\mathrm{mach}} \lVert D \rVert_F$ | Exact |
| **REQ-008** | Unitary Action Invariance | $\le c_4 M^2 \epsilon_{\mathrm{mach}} \lvert S(D) \rvert$ | $c_4 \approx 20$ |
| **REQ-008b** | Weyl Eigenvalue Bound | $\le c_5 M \epsilon_{\mathrm{mach}} \lVert D \rVert_2$ | $c_5 \approx 5$ |
| **REQ-009** | Central Finite Difference | $\le c_6 \epsilon_{\mathrm{mach}}^{2/3} \lVert \nabla S \rVert$ | $c_6 \approx 50$ |

### 3.4 Dual-Path Verification & Physical Invariants in Scientific Computing

#### 1. Verification Without Analytical Solutions
Computational physics software often lacks analytical solutions for general parameter regimes.
Interactive multi-matrix models with $g_4 > 0$ and dimension $N > 1$ have no closed-form solutions.
Standard unit tests cannot confirm whether an MCMC simulation samples the true physical distribution.
To solve this problem, leading scientific codes use cross-path verification and physical scaling laws.
Two independent algorithms compute the same physical quantity using different mathematical formulations.
Exact agreement between these independent pathways verifies algorithmic and physical correctness.

#### 2. Internal Energy Invariant: Fast Local Updates vs Global Traces
In statistical physics, the action $S(D)$ serves as the dimensionless potential energy of the noncommutative geometry.
Metropolis-Hastings MCMC accepts or rejects candidate configurations based on the internal energy difference:


$$
\Delta S = S(D_f) - S(D_i)
$$


The transition probability follows the Boltzmann factor:


$$
P_{\mathrm{accept}} = \min\left(1, \, \mathrm{e}^{-\Delta S}\right)
$$


RFL computes this quantity through two independent pathways:
* **Global Evaluation Path:** Reassembles the complete Dirac operator $D \in \mathbb{C}^{M \times M}$ and computes explicit matrix powers and traces in $\mathcal{O}(M^3)$ operations.
* **Local Incremental Path (`delta24`):** Exploits single-element matrix variations. It evaluates $\Delta S$ in $\mathcal{O}(N)$ operations using precomputed Clifford trace tensors ($\Omega$ table).

*Physical Failure Mode:*
If a sign error or indexing flaw corrupts the precomputed $\Omega$ tensor, the program does not crash.
All matrix configurations remain strictly Hermitian.
However, the Metropolis filter evaluates an incorrect energy difference.
The Markov chain drifts away from the true Boltzmann distribution $\mathrm{e}^{-S(D)}$ and samples an unphysical ensemble.
Verifying REQ-001 guarantees that the local update exactly mirrors the global internal energy change.

*Precedents in Major Scientific Software:*
* **Lattice QCD (Grid, Chroma, USQCD, MILC):**
  Lattice gauge updates evaluate the local action change using products of neighbouring link variables called staples.
  Lattice QCD suites maintain regression tests comparing local single-link staple updates against full 4D lattice action recalculations.
  They document this check as staple consistency or local-versus-global action invariance.
* **Molecular Dynamics (GROMACS, LAMMPS):**
  In microcanonical ($NVE$) ensembles, total energy $E = E_{\mathrm{kin}} + E_{\mathrm{pot}}$ is a conserved Hamiltonian invariant.
  GROMACS validates integrators using the `physical_validation` framework.
  The framework tracks energy drift against the shadow Hamiltonian to confirm physical correctness.

#### 3. Derivative Invariants: Generalised Forces & Symplectic Consistency
In physical dynamics, the negative gradient $-\nabla S$ represents the generalised force driving degrees of freedom.
Hamiltonian Monte Carlo (HMC) and Langevin algorithms integrate Hamilton's equations of motion:


$$
\dot{p}_k = -\frac{\partial S}{\partial M_k}
$$


RFL computes and verifies action gradients through two independent pathways:
* **Analytical Path:** Calculates symbolic matrix variations derived from trace cyclicity and noncommutative differential calculus.
* **Numerical Finite-Difference Path:** Evaluates directional derivatives using central finite difference stencils:


$$
\left( \nabla S \right)_{IJ} \approx \frac{S(M_{IJ} + h) - S(M_{IJ} - h)}{2h}
$$


*Physical Failure Mode:*
An error in analytical forces violates Liouville's theorem of phase space volume conservation.
In HMC, trajectories fail energy conservation ($\Delta H \gg 0$), collapsing Metropolis acceptance rates.
In Langevin simulations without acceptance steps, erroneous forces silently bias stationary expectation values.
Verifying REQ-009 guarantees symplectic consistency and correct force fields before executing dynamical sampling.

*Precedents in Major Scientific Software:*
* **Stan (Bayesian Modelling & HMC):**
  Stan provides an automated gradient verification diagnostic (`diagnose test=gradient`).
  Stan evaluates algorithmic automatic differentiation against numerical finite differences across random parameter vectors.
  Stan flags parameters exceeding a relative tolerance of $10^{-6}$ as algorithmic implementation bugs.
* **Lattice QCD (USQCD, Grid):**
  In dynamical pseudofermion simulations, the fermion force requires differentiating the inverted Dirac operator.
  USQCD suites maintain dedicated fermion force tests (`Test_fermion_force`).
  These tests integrate analytical forces along momentum trajectories and compare results against numerical action shifts.
* **Quantum Chemistry (PySCF, Psi4):**
  PySCF and Psi4 verify analytical Hellmann-Feynman gradients against numerical finite-difference energy derivatives.
  Every gradient module includes regression tests asserting agreement within $10^{-6}$ Hartree per Bohr.

#### 4. Scientific Documentation & Recording Standards
Leading scientific software projects document verification invariants through structured technical standards:
* **Unit-Test-First Coverage:** Whenever an invariant can be tested deterministically on a single state, projects implement it as a fast unit test.
* **Dual-Path Verification Tests:** Projects maintain automated regression tests executing both computational pathways in continuous integration.
* **Explicit Precision Budgets:** Documentation justifies numerical tolerances using machine precision and condition numbers instead of arbitrary thresholds.
* **Diagnostic Verification Tools:** Software exposes diagnostic commands (such as Stan's `test_grad` and GROMACS's `gmx check`) allowing users to verify algorithmic consistency.
* **Theory-to-Code Traceability:** Technical documentation links code routines directly to underlying theoretical equations and literature citations.

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
* **Tier 1 (Smoke / Deterministic Gate, < 2 minutes):** Executes on every commit and PR. Enforces the unit-test-first policy. Tests all deterministic algebraic identities, Hermiticity, gauge invariance, and finite differences in fast unit tests.
* **Tier 2 (Physics Integration Gate, ~3–5 minutes):** Executes on PR merge to `main`. Uses fixed seeds to test short MCMC chains, loose-tolerance Gaussian moments ($4\sigma$), and state-flux balance.
* **Tier 3 (Scientific Validation Suite, Nightly / Release Gate, ~30–60 minutes):** Executes high-statistics runs, multi-chain split $\hat{R}$, Riemann-Hilbert curve fits, and Simulation-Based Calibration.
* **Flakiness Control Protocol:** Stochastic tests apply the Holm-Bonferroni correction to prevent false discovery inflation. If a stochastic test fails with $p < \alpha$, CI triggers an automated second run with an independent seed and double chain length before failing.

---

## 5. Target Architecture & Component Contracts

### 5.1 C++ Core Verification (`src/core/tests/`)
C++ tests focus on deterministic invariants and algebraic correctness:
* `tDelta.cpp`: Verifies incremental action update `delta24` against full action difference using dimension-scaled tolerance:
  ```cpp
  TEST(DeltaTests, IncrementalActionMatchesFullDifference) {
    const DiracOperator dirac(1, 3, 6);
    const Action action(-1.0, 1.0);
    const double Si = action.calculateS(dirac);
    const double dS = dirac.delta24(action, 0, 1, 2, {0.1, 0.0});
    // Apply update and recalculate Sf
    const double Sf = action.calculateS(dirac);
    const double M = dirac.getMatrixDimension() * dirac.getGammaDimension();
    const double tol = 10.0 * M * std::numeric_limits<double>::epsilon() * (std::abs(Sf) + std::abs(Si)) + 1e-12;
    EXPECT_NEAR(Sf - Si, dS, tol);
  }
  ```
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
└── test_scaling_invariants.py
```

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Test Description | Target Invariant | Target Tier |
| :--- | :--- | :--- | :--- |
| **Algebraic Invariants** | `ctest --test-dir build -R DeltaTests` | Verifies $\Delta S$ matches full action differences (REQ-001). | Tier 1 (PR Gate) |
| **Spectral Axioms** | `pytest tests/physics/test_axiomatic_invariants.py` | Verifies Hermiticity, KO table, and chirality (REQ-002–004). | Tier 1 (PR Gate) |
| **Gauge Invariance** | `pytest tests/physics/test_gauge_invariance.py` | Verifies $S(UDU^\dagger) = S(D)$ and spectral invariance (REQ-008). | Tier 1 (PR Gate) |
| **Derivative Consistency** | `pytest tests/physics/test_derivative_consistency.py` | Verifies variations match central finite differences (REQ-009). | Tier 1 (PR Gate) |
| **Gaussian Limit** | `pytest tests/physics/test_gaussian_limit.py` | Verifies convergence to Wigner convolution law (REQ-005). | Tier 2 (Merge Gate) |
| **Detailed Balance** | `pytest tests/physics/test_detailed_balance.py` | Verifies microscopic reversibility and state-flux balance (REQ-006). | Tier 2 (Merge Gate) |
| **Physical Scaling Laws** | `pytest tests/physics/test_scaling_invariants.py` | Verifies parameter rescaling and coupling monotonicity (REQ-011). | Tier 2 (Merge Gate) |
| **1-Matrix Benchmark** | `pytest tests/physics/test_riemann_hilbert_1matrix.py` | Verifies agreement with Riemann-Hilbert analytical solutions. | Tier 3 (Nightly) |
| **MCMC Calibration (SBC)** | `pytest tests/physics/test_gaussian_sbc.py` | Verifies uniform posterior ranks in Gaussian ensembles (REQ-010). | Tier 3 (Nightly) |
| **Chain Convergence** | `pytest tests/physics/test_chain_convergence.py` | Verifies rank-normalized split $\hat{R} < 1.05$ and $\mathrm{ESS} \ge 400$ (REQ-007). | Tier 3 (Nightly) |

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
  1. Implement `tests/physics/test_riemann_hilbert_1matrix.py` verifying signature $(0, 1)$ spectral densities against exact analytical curves.
  2. Implement `tests/physics/test_detailed_balance.py` testing microscopic reversibility and coarse-grained state-flux balance.
  3. Implement `tests/physics/test_gaussian_sbc.py` validating the MCMC sampler with Simulation-Based Calibration.
  4. Implement rank-normalized folded split $\hat{R}$ and $\mathrm{ESS}_{\mathrm{bulk}}$ diagnostics in Python analysis utilities.
  5. Implement `tests/physics/test_scaling_invariants.py` to verify parameter rescaling and coupling monotonicity.
  6. Create benchmark scripts to recreate published matrix model scaling curves and archive golden reference datasets.
  7. Configure automated Zenodo DOI archiving upon GitHub release tags.
