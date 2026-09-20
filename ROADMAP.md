# RFL Development Roadmap

> [!IMPORTANT]
> **Active Milestone:** `v0.3.0` (Core Modernisation & Scientific Verification)  
> **Current Focus:** `Phase 0: Pre-flight & Baseline Assurance`  
> **Milestone Tracker:** [Issue #71](https://github.com/pauldruce/RFL/issues/71) | [GitHub Project Board](https://github.com/users/pauldruce/projects/1)

This document defines the release roadmap, execution order, and component dependencies for the Random Fuzzy Library (RFL).

---

## 1. Release Progression Overview

| Release | Primary Focus | Status |
| :--- | :--- | :--- |
| **`v0.1.0`** | Python wheel packaging, `pybind11` bindings, Carma zero-copy NumPy interoperability, ASD-STE100 standardisation | ✅ Shipped |
| **`v0.2.0`** | Dependency decoupling with CMake `FetchContent`, tiered CI pipeline, native Windows MSVC support | ✅ Shipped |
| **`v0.3.0`** | Value-type Dirac operator, non-blocking simulation stepper API, CMake target exports, and scientific verification suite | 🔄 Active |
| **`v0.4.0`** | Community package managers (Homebrew, Conda-Forge, Conan, vcpkg) and asymptotic Riemann-Hilbert benchmarks | 💡 Planned |

---

## 2. Milestone v0.3.0 Phased Delivery Plan

The `v0.3.0` milestone executes across five sequential phases.
Each phase builds on the verified invariants of previous phases to prevent regressions.

```
Phase 0: Pre-flight & Baseline Assurance
  ├── #62 (Windows vcpkg DLL copying race)
  ├── #56 (Codify C++ & Python contract baseline tests)
  ├── #67 (Replace GSL with StdRng, BSD-3-Clause licence, CITATION.cff)
  └── #49 (Fix MSVC /W4 conversion warnings & discrete sampling)

Phase 1: Core Mathematical Foundation
  ├── #46 & #12 (DiracOperator value type, lazy omega table, strict const)
  ├── #14 (Extract BarrettGlaserAction with delta24 & gradients)
  └── EP-4 Task 1 (Port tDelta.cpp scale-aware unit tests)

Phase 2: Sampler Decoupling & Stepper Engine
  └── #59 & #13 (Raw pointer elimination, non-blocking Stepper API, dual-averaging)

Phase 3: Scientific Verification Suite
  └── #69 (Axioms, gauge invariance, finite differences, Gaussian limit, CI tiers)

Phase 4: Build System, Packaging & Distribution
  ├── #60 (Adopt canonical include/rfl/ layout)
  ├── #70 (Modular CMake targets and Precompiled Headers)
  └── #3 (CMake install targets, RFLConfig.cmake, CPack archives)

Phase 5: Documentation, Profiling & Release Qualification
  ├── #57 (Doxygen API docs & GitHub Pages CI)
  ├── #58 (Armadillo / OpenMP multi-threading investigation)
  └── Release Qualification & Tagging
```

---

### Detailed Phase Breakdown

| Phase | Component | Issues & Deliverables | Dependencies | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 0** | **Pre-flight & Quality Assurance** | • [#62](https://github.com/pauldruce/RFL/issues/62): Fix Windows parallel DLL file collisions<br/>• [#56](https://github.com/pauldruce/RFL/issues/56): Compile-time C++ contract and Python API surface tests<br/>• [#67](https://github.com/pauldruce/RFL/issues/67): Implement `StdRng`, remove GSL, add BSD-3-Clause `LICENSE` & `CITATION.cff`<br/>• [#49](https://github.com/pauldruce/RFL/issues/49): Address MSVC conversion warnings and uniform discrete sampling | None | 🔄 In Progress |
| **Phase 1** | **Core Mathematical Foundation** | • [#46](https://github.com/pauldruce/RFL/issues/46): Optimise `DiracOperator` omega table initialisation<br/>• [#12](https://github.com/pauldruce/RFL/issues/12): Refactor `DiracOperator` into regular value type with strict `const` correctness<br/>• [#14](https://github.com/pauldruce/RFL/issues/14): Extract analytic trace variations into `BarrettGlaserAction`<br/>• Port `src/core/tests/tDelta.cpp` unit tests | Phase 0 | ⏳ Scheduled |
| **Phase 2** | **Sampler Decoupling & Stepper API** | • [#59](https://github.com/pauldruce/RFL/issues/59): Eliminate raw pointer members in `Metropolis` and `Hamiltonian`<br/>• [#13](https://github.com/pauldruce/RFL/issues/13): Implement `MetropolisSampler` stepper, `ISimulationObserver`, and dual-averaging step calibration | Phase 1 | ⏳ Scheduled |
| **Phase 3** | **Scientific Verification Suite** | • [#69](https://github.com/pauldruce/RFL/issues/69): Implement EP-4 Phase 1 test suite (Hermiticity, 8-fold KO table, chirality, Haar gauge invariance, finite differences, Gaussian limit moments, CI tiers) | Phase 2 | ⏳ Scheduled |
| **Phase 4** | **Build & Package Distribution** | • [#60](https://github.com/pauldruce/RFL/issues/60): Migrate headers to canonical `include/rfl/` structure<br/>• [#70](https://github.com/pauldruce/RFL/issues/70): Modularise CMake targets and add Precompiled Headers<br/>• [#3](https://github.com/pauldruce/RFL/issues/3): CMake `install()` rules, `RFLConfig.cmake`, and CPack binary packaging | Phase 1 | ⏳ Scheduled |
| **Phase 5** | **Documentation & Polish** | • [#57](https://github.com/pauldruce/RFL/issues/57): Automate Doxygen reference and GitHub Pages deployment<br/>• [#58](https://github.com/pauldruce/RFL/issues/58): Investigate Armadillo and OpenMP multi-threading concurrency<br/>• Release qualification and PyPI distribution | Phase 3, 4 | ⏳ Scheduled |

---

## 3. Guiding Enhancement Proposals

Every architectural modernisation task traces directly to an approved Enhancement Proposal in `docs/eps/`:

* **[EP-1: Core Architecture Modernisation](docs/eps/ep-1-core-architecture-modernisation.md):** Governs `DiracOperator` value semantics, the non-blocking Stepper API, and `BarrettGlaserAction` decoupling.
* **[EP-2: Multi-Platform Package and Binary Distribution](docs/eps/ep-2-package-and-binary-distribution.md):** Governs CMake `install()` targets, `RFLConfig.cmake` package files, and CPack release archives.
* **[EP-3: Resource-Efficient CI Pipeline](docs/eps/ep-3-resource-efficient-ci-pipeline.md):** Governs modular compilation targets, build caching, and Precompiled Headers.
* **[EP-4: Computational Physics Verification Suite](docs/eps/ep-4-computational-physics-verification-suite.md):** Governs physical validation standards, geometric axioms, gauge symmetries, and limiting spectral distributions.

---

## 4. Contributing & Working on a Phase

1. Always check the **Milestone Tracker** ([#71](https://github.com/pauldruce/RFL/issues/71)) and the [Project Board](https://github.com/users/pauldruce/projects/1) before starting work.
2. Ensure upstream dependencies for the target phase are completed and merged to `main`.
3. Create feature branches using conventional naming: `feat/issue-<N>-<short-name>` or `fix/issue-<N>-<short-name>`.
4. Reference the issue in PR descriptions using `Closes #<N>`.
