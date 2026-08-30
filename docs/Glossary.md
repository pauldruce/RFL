# RFL Controlled Vocabulary & Glossary

This document establishes the official **ASD-STE100 Controlled Vocabulary** for RFL. 

### The Core Principle: One Word, One Meaning
> *Each approved term has exactly one defined meaning. Do not use synonyms for the same concept. Do not use the same word for different concepts.*

All documentation, code docstrings, comments, Enhancement Proposals (EPs), and research notes must use the approved terms defined in this glossary.

---

## 1. Physical & Mathematical Terms

| Approved Term | Definition | Unapproved Synonyms (Avoid) |
| :--- | :--- | :--- |
| **Dirac Operator** | The matrix operator $D = \sum \gamma^i \otimes H_i + \sum \gamma^{p+j} \otimes L_j$ that defines the geometry. | *Dirac matrix (when referring to the operator object)* |
| **Clifford Module** | The real or complex representation of the Clifford algebra $\text{Cl}(p,q)$ holding gamma matrices $\gamma^a$. | *Clifford algebra engine, Gamma set* |
| **Gamma Matrix** | A generator matrix $\gamma^a$ of the Clifford algebra satisfying $\{\gamma^a, \gamma^b\} = 2\eta^{ab}\mathbb{I}$. | *Dirac gamma, Clifford matrix* |
| **Spectral Triple** | The mathematical triplet $(A, \mathcal{H}, D)$ that defines a noncommutative geometry. | *NCG triplet, Connes triple* |
| **Barrett-Glaser Action** | The spectral action functional $S(D) = g_2 \text{Tr}(D^2) + g_4 \text{Tr}(D^4)$. | *BG functional, Spectral potential* |
| **Signature** | The integer pair $(p, q)$ representing $p$ Hermitian and $q$ anti-Hermitian matrix generators. | *Clifford type, metric signature* |
| **Matrix Dimension ($N$)** | The size $N \times N$ of the internal matrices $H_i, L_j$. | *Matrix size, cutoff size* |
| **Chirality Operator ($\Gamma$)** | The grading operator satisfying $\Gamma^2 = \mathbb{I}$ and $\{\Gamma, D\} = 0$ for even geometries. | *Grading operator, gamma five* |
| **Reality Operator ($J$)** | The antilinear isometry representing charge conjugation in real spectral triples. | *Charge conjugation operator* |

---

## 2. Computational & Sampling Terms

| Approved Term | Definition | Unapproved Synonyms (Avoid) |
| :--- | :--- | :--- |
| **Step** | A single coordinate proposal and Metropolis accept/reject test on one matrix element. | *Iteration, move, jump* |
| **Sweep** | A complete pass of steps across all matrix degrees of freedom ($N_{\text{mat}} \times N^2$). | *Epoch, generation, cycle* |
| **Thermalisation** | The initial sampling phase to reach statistical equilibrium before recording data. | *Burn-in, warmup, equilibration* |
| **Dual-Averaging** | The adaptive algorithm that tunes the proposal step scale to reach a target acceptance rate. | *NUTS tuning, auto-scaling* |
| **Acceptance Rate** | The ratio of accepted proposals to total proposals in a sweep ($\in [0, 1]$). | *Acceptance ratio, hit rate* |
| **Observable** | A physical quantity calculated from the Dirac operator (e.g. $\text{Tr}(D^2)$, $\text{Tr}(D^4)$, eigenvalues). | *Measurement, metric, output* |
| **Observer** | A software object that listens to step/sweep events and records observables. | *Sink, listener, telemetry logger* |
| **Autocorrelation Time ($\tau_{\text{int}}$)** | The integrated statistical correlation time between successive Markov samples. | *Correlation length, memory time* |
| **Eigenvalue Spectrum** | The set of real eigenvalues $\{\lambda_i\}$ of the assembled Dirac operator $D$. | *Eigen-spectrum, Dirac energies* |

---

## 3. Approved Technical Verbs (ASD-STE100)

| Approved Verb | Approved Meaning & Usage | Unapproved Words (Avoid) |
| :--- | :--- | :--- |
| **Calculate** | To determine a value using exact mathematical equations or analytical formulas. | *Compute (when analytical), evaluate (when vague)* |
| **Compute** | To determine a value using iterative numerical algorithms or matrix solvers. | *Calculate (when iterative), figure out* |
| **Initialise** | To set up starting memory, configurations, or objects before use. | *Instantiate, boot up, prepare* |
| **Randomise** | To populate matrix elements with random values drawn from an RNG engine. | *Scramble, randomize (US spelling)* |
| **Verify** | To test and prove that a state satisfies a mathematical invariant or requirement. | *Check, ensure, validate (when testing invariants)* |
| **Sample** | To draw Markov chain states from the probability distribution $e^{-S(D)}$. | *Simulate (when referring to the MCMC step)* |
| **Run** | To start and execute a simulation process or test suite. | *Carry out, perform, execute* |
| **Stop** | To bring an iteration or simulation process to an immediate end. | *Terminate, abort, kill, halt* |

---

## 4. Controlled Conjunctions & Prepositions

| Approved Word | Usage Rule | Forbidden Alternatives |
| :--- | :--- | :--- |
| **`because`** | Use only to give a logical reason. | *as, since, for the reason that* |
| **`when`** | Use only to specify a condition or simultaneous event. | *as, at the time that* |
| **`after`** | Use only to specify a sequence in time. | *since, following on from* |
| **`to`** | Use to indicate purpose. | *in order to, so as to* |
| **`must`** | Use only for mandatory requirements. | *shall, should, ought to* |
| **`can`** | Use only for capability or possibility. | *may (when meaning capability)* |

---

## 5. British vs. US Spelling Dictionary

RFL strictly enforces **British English** spelling:

| British English (Approved) | US English (Forbidden) | Example Context |
| :--- | :--- | :--- |
| **Initialise** | *Initialize* | `dirac.initialise()` |
| **Randomise** | *Randomize* | `dirac.randomise(rng)` |
| **Optimise / Optimisation** | *Optimize / Optimization* | Compiler optimisation |
| **Diagonalisation** | *Diagonalization* | Exact matrix diagonalisation |
| **Behaviour** | *Behavior* | Scaling behaviour |
| **Modelling** | *Modeling* | Matrix modelling |
| **Centre** | *Center* | Spectral centre |
| **Programme** | *Program* | Research programme *(use 'program' only for computer code)* |
| **Fibre** | *Fiber* | Clifford fibre bundle |
| **Gauge** | *Gage* | Gauge field fluctuations |
| **Analogue** | *Analog* | Continuum analogue |
