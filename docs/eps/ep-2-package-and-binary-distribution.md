# EP-2: Multi-Platform Package and Binary Distribution

* **Title:** Multi-Platform Package and Binary Distribution
* **Author:** Paul Druce
* **Status:** Accepted (In Progress)
* **Target Versions:** RFL v0.1.0 (Phase 1), v0.3.0 (Phase 2), Backlog (Phase 3)
* **Date:** 2026-08-30

---

## 1. Multi-Phase Implementation Tracker

This proposal spans multiple releases.
The table below tracks the status of each implementation phase:

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Directory layout, PyPI binary wheels, and release CI/CD | `v0.1.0` | [PR #30](https://github.com/pauldruce/RFL/pull/30) (Closes [#18](https://github.com/pauldruce/RFL/issues/18)) | ✅ Implemented |
| **Phase 2** | CMake `install()` targets, `RFLConfig.cmake`, and CPack archives | `v0.3.0` | [#3](https://github.com/pauldruce/RFL/issues/3) | ⏳ Scheduled |
| **Phase 3** | Package manager distribution (Homebrew Tap, Conda-Forge, Conan) | Future | [#29](https://github.com/pauldruce/RFL/issues/29) | 💡 Planned |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Research Context
`RFL` enables researchers to simulate finite noncommutative geometries and random fuzzy spaces using Markov Chain Monte Carlo methods.
Researchers study spectral triples $(A, \mathcal{H}, D)$ and the Barrett-Glaser action across diverse computing environments:
1. Python data science workflows (Jupyter notebooks, NumPy, SciPy, Matplotlib).
2. High-performance C++ batch simulations.
3. Automated cluster runs on Linux and macOS workstations.

Previously, consuming RFL required compiling C++ source code locally.
Consumers had to install CMake, C++17 compilers, Armadillo, LAPACK, BLAS, and GSL manually.
This requirement created friction for physicists and data scientists who require a quick `pip install pyrfl` workflow.
In addition, external C++ simulation programs had no standard CMake target mechanism to link `rfl`.
The codebase also suffered from historical directory naming debt (`new_RFL` and `old_RFL`), which created confusing target names.

### 2.2 Goals
* **Clean Target & Directory Standardisation:** Promote `new_RFL` to `core/` with target `rfl_core` (aliases: `RFL::core`, `RFL::rfl`), and archive `old_RFL` to `legacy/` (`rfl_legacy`).
* **Multi-Platform Binary Wheels:** Provide precompiled, standalone Python wheels on PyPI for Linux (`manylinux_2_28`) and macOS (`x86_64`, `arm64`) across Python 3.9–3.13 (`pip install pyrfl`).
* **Hermetic Dynamic Library Vendoring:** Automatically vendor and bundle shared C++ runtime dependencies (`openblas`, `gsl`, `armadillo`) inside wheels via `auditwheel` and `delocate`.
* **Zero-Secret CI/CD PyPI Publishing:** Implement automated PyPI deployment on version tag push using OpenID Connect (OIDC) Trusted Publishing.
* **CMake In-Tree FetchContent Support:** Ensure external C++ codebases can integrate RFL seamlessly via standard `FetchContent`.
* **Foundation for Package Managers:** Lay the ground for future CMake `install()` targets, CPack release archives, and community package managers (Homebrew, Conda-Forge).

### 2.3 Non-Goals
* Supporting legacy Python versions (< 3.9) or non-standard interpreters (PyPy).
* Supporting 32-bit operating systems or musl libc distributions.
* Distributing closed-source binary blobs without public build recipes.

---

## 3. Research Workflows & Scientific Requirements

### 3.1 Core Research Scenarios
1. **Scenario 1 (Interactive Spectral Analysis in Python):** A physicist installs RFL via `pip install pyrfl` in an isolated virtual environment and computes Dirac eigenvalue spectra in Jupyter (`import rfl`).
2. **Scenario 2 (Custom C++ Simulation Pipeline):** A researcher includes RFL into an external C++ MCMC solver using `FetchContent_Declare(RFL ...)` and links `RFL::core` without manually configuring header search paths.
3. **Scenario 3 (Unified Scientific Environment):** A research lab installs `librfl` and `rfl` via Conda into reproducible simulation environments.

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Physical & Technical Invariant |
| :--- | :--- | :--- |
| **REQ-PKG-01** | **Directory & Target Standardisation** | Restructure `new_RFL/` → `core/` (target `rfl_core`, alias `RFL::core`) and `old_RFL/` → `legacy/` (target `rfl_legacy`). |
| **REQ-PKG-02** | **Multi-Platform Binary Wheels** | Build standalone wheels for Linux (`manylinux_2_28_x86_64`) and macOS (`macosx_14_0_arm64`, `macosx_15_0_x86_64`) covering Python 3.9–3.13. |
| **REQ-PKG-03** | **Vendored Shared Libraries** | Bundle dynamic dependencies (`openblas`, `gsl`, `armadillo`) so wheels execute on clean systems. |
| **REQ-PKG-04** | **Pre-Publication Test Gate** | Execute Python test suite (`pytest`) inside clean wheel environments before publication. |
| **REQ-PKG-05** | **Secure PyPI Publishing** | Use OpenID Connect (OIDC) Trusted Publishing to prevent static secret exposure. |
| **REQ-PKG-06** | **CMake FetchContent Support** | Export namespaced alias `RFL::core` and `RFL::rfl` in top-level CMake configuration. |
| **REQ-PKG-07** | **Standard CMake Installation** | Define CMake `install()` targets for public headers, compiled libraries, and CMake config packages. |
| **REQ-PKG-08** | **CPack Binary Packaging** | Generate `.tar.gz` release archives containing headers, static libraries, and package configuration files. |
| **REQ-PKG-09** | **Ecosystem Package Distribution** | Provide recipes for Homebrew Tap, Conda-Forge, and evaluate C++ registries (Conan, vcpkg). |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: Directory Layout & Target Renaming

| Criteria | Option A: Standardised Layout (`core/` & `legacy/`) (Selected) | Option B: Retain `new_RFL` / `old_RFL` |
| :--- | :--- | :--- |
| **Public API Clarity** | High (`RFL::core`, `RFL::rfl`) | Poor (`RFL::new_RFL` conveys prototype state) |
| **Historical Separation** | Clear (`legacy/` explicitly marked) | Ambiguous |
| **Packaging Cleanliness** | Standard CMake conventions | Non-standard paths |
| **Decision** | **Selected (Option A)** | Rejected |

*Rationale:* Promoting modern code to `src/RFL/core` and archiving historical code to `src/RFL/legacy` establishes an intuitive and lasting naming standard before the `v0.1.0` release.

---

### 4.2 ADR-2: Python Wheel Build Automation Framework

| Criteria | Option A: `pypa/cibuildwheel` (Selected) | Option B: Custom Matrix Shell Scripts | Option C: Source Distribution (`sdist`) Only |
| :--- | :--- | :--- | :--- |
| **Standardisation** | High (PyPA standard for scientific Python) | Low (Bespoke maintenance) | High (Standard `build --sdist`) |
| **Library Vendoring** | Automated (`auditwheel`, `delocate`) | Manual configuration | None (Fails on systems without C++ stack) |
| **Isolated Testing** | Built-in per Python version | Custom virtual environments | Relies on consumer environment |
| **User Experience** | Instant `pip install pyrfl` | Instant `pip install pyrfl` | High friction for non-C++ users |
| **Decision** | **Selected (Option A)** | Rejected | Rejected |

*Rationale:* `cibuildwheel` encapsulates `manylinux` container management, executes tests inside clean wheel environments, and automates dynamic library vendoring.

---

### 4.3 ADR-3: C++ Dependency Management & Library Vendoring

| Criteria | Option 1: Package Manager Provisioning + Vendoring (Selected) | Option 2: Monolithic Static Build from Source |
| :--- | :--- | :--- |
| **CI Build Time** | Fast (1–3 minutes using `dnf`/`brew`) | Slow (15–30 minutes compiling OpenBLAS/GSL) |
| **Binary Portability** | High (`auditwheel` and `delocate` bundle shared libraries) | High (Static linking) |
| **Maintenance** | Low (System package updates) | High (Custom CMake build recipes) |
| **Decision** | **Selected (Option 1)** | Rejected |

*Rationale:* Installing `openblas`, `lapack`, `gsl`, and `armadillo` via `dnf` on Linux and `brew` on macOS provides fast, reliable CI builds.
`auditwheel` and `delocate` automatically discover and vendor dynamic dependencies into the wheel.
On macOS, `MACOSX_DEPLOYMENT_TARGET` is explicitly set to `14.0` on Apple Silicon (`arm64`) and `15.0` on Intel (`x86_64`) to match Homebrew binary SDK baselines.

---

### 4.4 ADR-4: PyPI Authentication & Supply Chain Security

| Criteria | Option A: PyPI Trusted Publishing / OIDC (Selected) | Option B: Repository Secret API Token |
| :--- | :--- | :--- |
| **Security Posture** | High (Short-lived cryptographic tokens) | Moderate (Long-lived secret stored in GitHub) |
| **Secret Management** | Zero secrets required in repository settings | Requires updating static tokens |
| **Industry Standard** | PyPA recommended standard | Legacy method |
| **Decision** | **Selected (Option A)** | Rejected |

*Rationale:* Trusted Publishing authenticates directly between GitHub Actions and PyPI using OpenID Connect, eliminating credential leakage risks.

---

### 4.5 ADR-5: C++ Consumption & Distribution Strategy

| Criteria | Option A: CMake `FetchContent` (Phase 1) | Option B: CPack Binary Tarballs (Phase 2) | Option C: Package Registries (Phase 3) |
| :--- | :--- | :--- | :--- |
| **Target Audience** | Active C++ researchers building from source | C++ researchers wanting prebuilt libraries | macOS and Conda data scientists |
| **ABI Compatibility** | Perfect (Builds with consumer compiler) | Restricted to compatible toolchains | Managed by package manager |
| **Infrastructure** | Zero extra infrastructure | GitHub Releases asset attachment | Feedstock / Tap repositories |
| **Decision** | **Selected (Phase 1)** | **Selected (Phase 2)** | **Selected (Phase 3)** |

*Rationale:* Combining all three distribution channels across distinct phases provides immediate integration while establishing a robust long-term packaging roadmap.

---

### 4.6 ADR-6: Package Manager Ecosystem Trade-offs

| Package Manager | Target Ecosystem | Pros | Cons & Considerations |
| :--- | :--- | :--- | :--- |
| **Homebrew** | macOS Researchers | Seamless `brew install rfl` resolving dependencies. | macOS focus; Linuxbrew adoption is low. |
| **Conda-Forge** | Scientific Python & C++ | Packages both `librfl` and `rfl` with pinned BLAS/GSL. | Requires maintaining a separate feedstock repo. |
| **vcpkg / Conan** | Enterprise & C++ Developers | Integrated C++ dependency graph management. | Higher overhead; less common in mathematical physics. |

---

## 5. Target Architecture & Workflow Design

### 5.1 The 3-Phase Roadmap

```mermaid
flowchart TD
    P1["Phase 1 (v0.1.0)
    • Core / Legacy layout
    • cibuildwheel (Linux & macOS)
    • PyPI OIDC Trusted Publishing
    • CMake FetchContent"]
    P2["Phase 2 (v0.3.0)
    • CMake install() targets
    • Export RFLConfig.cmake
    • CPack binary .tar.gz archives"]
    P3["Phase 3 (Backlog)
    • Homebrew Tap formula
    • Conda-Forge feedstock
    • Evaluate Conan & vcpkg"]

    P1 ==> P2 ==> P3
```

| Phase | Target Version | Scope & Key Deliverables | Ecosystem |
| :--- | :--- | :--- | :--- |
| **Phase 1** | `v0.1.0` | • Decompose `core/` and `legacy/`<br/>• `cibuildwheel` across Linux and macOS<br/>• PyPI OIDC Trusted Publishing<br/>• CMake `FetchContent` support (`RFL::core`) | Python & CMake FetchContent |
| **Phase 2** | `v0.3.0` | • CMake `install()` rules<br/>• Export `RFLConfig.cmake` package files<br/>• CPack binary `.tar.gz` distributions | System C++ Installations |
| **Phase 3** | Backlog | • Homebrew Tap for macOS (`brew install pauldruce/rfl/rfl`)<br/>• Conda-Forge feedstock for scientific Python<br/>• Conan / vcpkg packaging | Package Managers |

### 5.2 Release Pipeline Workflow (`.github/workflows/release.yml`)

```mermaid
flowchart TD
    Trigger(["Trigger: Release Published or Workflow Dispatch"])

    Trigger --> BuildWheels["cibuildwheel Matrix
    • Linux x86_64 (manylinux_2_28)
    • macOS x86_64 (macos-15-intel)
    • macOS arm64 (macos-14)"]
    Trigger --> BuildSdist["Source Distribution (sdist)
    • Canonical source tarball (.tar.gz)"]

    BuildWheels --> Stage["Staged Release Artifacts
    (dist/*.whl & dist/*.tar.gz)"]
    BuildSdist --> Stage

    Stage --> GH["GitHub Release
    (Upload Assets)"]
    Stage --> PyPI["PyPI Production
    (Trusted OIDC)"]
```

| Pipeline Stage | Job Name | Execution Environment | Deliverables & Verification |
| :--- | :--- | :--- | :--- |
| **1. Trigger** | `release` (published) or `workflow_dispatch` | GitHub Actions | Dispatches matrix jobs concurrently |
| **2. Binary Wheels** | `build_wheels` | • Linux x86_64 (`manylinux_2_28`)<br/>• macOS x86_64 (`macos-13`)<br/>• macOS arm64 (`macos-14`) | Compiles wheels via `cibuildwheel` and executes `pytest` |
| **3. Source Dist** | `build_sdist` | Ubuntu latest | Generates canonical `.tar.gz` sdist |
| **4. Asset Upload** | `upload_release_assets` | Ubuntu latest (on published release) | Attaches wheels and sdist to GitHub Release |
| **5. PyPI Publish** | `publish_pypi` | Ubuntu latest (OIDC Trusted Publishing) | Publishes packages to PyPI index |

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Test Description | Target Invariant |
| :--- | :--- | :--- |
| **C++ Core Unit Tests** | `ctest --test-dir build --output-on-failure` | Verifies `rfl_core` and `rfl_legacy` suites pass completely. |
| **Local Unit Tests** | `pytest src/RFL/python_bindings/tests/unit` | Verifies Python API and exception propagation. |
| **Local Integration Tests** | `pytest src/RFL/python_bindings/tests/integration` | Verifies NumPy zero-copy buffer safety and memory management. |
| **Please Build Verification** | `./pleasew test //src/RFL/python_bindings/tests:all` | Validates in-tree Please hermetic build with new `core/` path. |
| **cibuildwheel Test Phase** | `pytest {project}/src/RFL/python_bindings/tests` | Verifies installed wheel inside clean containers. |
| **Auditwheel Verification** | `auditwheel show <wheel_file>` | Verifies all dynamic libraries are vendored into `.libs`. |
| **Delocate Verification** | `delocate-listdeps <wheel_file>` | Verifies macOS dylibs are vendored and paths patched. |

---

## 7. Phased Delivery Plan

### Phase 1: Directory Restructuring, PyPI Wheels, GitHub Releases & CMake FetchContent
* **Target Version:** `v0.1.0`
* **GitHub Issue:** [Issue #18](https://github.com/pauldruce/RFL/issues/18)
* **Tasks:**
  1. Rename `src/RFL/new_RFL` → `src/RFL/core` and `src/RFL/old_RFL` → `src/RFL/legacy`.
  2. Update CMake target names (`rfl_core`, `RFL::core`, `RFL::rfl`, `rfl_legacy`) across CMakeLists.txt.
  3. Update Please build targets (`src/RFL/core/BUILD`, `src/RFL/legacy/BUILD`, `src/RFL/BUILD`).
  4. Configure `[tool.cibuildwheel]` in `src/RFL/pyproject.toml`.
  5. Create `.github/workflows/release.yml` with wheel matrix, sdist build, GitHub Release asset upload, and PyPI OIDC publishing.
  6. Update `README.md` and `docs/Consumption_Guide.md` with `pip` and `FetchContent` instructions.

### Phase 2: CMake Installation Targets & CPack Binary Packaging
* **Target Version:** `v0.3.0`
* **GitHub Issue:** [Issue #3](https://github.com/pauldruce/RFL/issues/3)
* **Tasks:**
  1. Define CMake `install(TARGETS rfl_core EXPORT RFLTargets ...)` and header install rules.
  2. Generate and install `RFLConfig.cmake` and `RFLConfigVersion.cmake`.
  3. Configure CPack to package `.tar.gz` and `.zip` archives.
  4. Add CPack artifact generation step to `.github/workflows/release.yml`.

### Phase 3: Community Package Managers
* **Target Version:** Future / Backlog
* **GitHub Issue:** [Issue #29](https://github.com/pauldruce/RFL/issues/29)
* **Tasks:**
  1. Create Homebrew formula in `pauldruce/homebrew-rfl`.
  2. Create Conda-Forge feedstock for `librfl` and `rfl`.
  3. Complete research trade study on Conan and vcpkg registries.
