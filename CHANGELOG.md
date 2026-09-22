# Changelog

All notable changes to the Random Fuzzy Library (RFL) are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.3.0] - Unreleased

### Highlights & Breaking Changes
* **Breaking Change (Core GSL Removal & `StdRng`):** Replaced GNU Scientific Library (GSL) in `rfl_core` with standard C++ `StdRng` based on `std::mt19937_64`. GSL is now optional and isolated strictly to historical `rfl_legacy` builds ([#67](https://github.com/pauldruce/RFL/issues/67)).
* **PyPI Wheel Decoupling:** Stripped `libgsl` dependency and copyleft constraints from precompiled binary Python wheels (`cibuildwheel`) and CI runner workflows ([#67](https://github.com/pauldruce/RFL/issues/67)).
* **Formal BSD-3-Clause Licencing:** Added repository root `LICENSE` (BSD-3-Clause) and declared PEP 639 `license = "BSD-3-Clause"` in `pyproject.toml` ([#67](https://github.com/pauldruce/RFL/issues/67)).
* **Academic Citation Metadata:** Added root `CITATION.cff` conforming to Citation File Format 1.2.0 for 1-click citation export on GitHub, Zenodo, and Google Scholar ([#67](https://github.com/pauldruce/RFL/issues/67)).
* **API Contract & Surface Testing:** Enforced compile-time C++ contract tests and runtime Python API surface verification ([#74](https://github.com/pauldruce/RFL/pull/74)).
* **Codebase Cleanliness & Modern Tooling:** Fixed header guard collisions, removed dead commented legacy tests, and repaired verification scripts ([#75](https://github.com/pauldruce/RFL/pull/75)).

### 🚀 Features & Enhancements
* refactor(core): replace GSL with standard C++ random engine and adopt BSD-3-Clause licence ([#67](https://github.com/pauldruce/RFL/issues/67)) by @pauldruce in [#77](https://github.com/pauldruce/RFL/pull/77)
* test(api): establish compile-time C++ contract tests and Python API surface assurance by @pauldruce in [#74](https://github.com/pauldruce/RFL/pull/74)

### 🐛 Bug Fixes
* fix(ci,windows): resolve concurrent vcpkg z-applocal DLL copying race condition by @pauldruce in [#73](https://github.com/pauldruce/RFL/pull/73)
* fix(ci): repair status badges in README and declare packaging classifiers by @pauldruce in [#65](https://github.com/pauldruce/RFL/pull/65)

### 🧰 Build & CI/CD Architecture
* chore(clean): fix header guard collisions, clean commented legacy tests, and repair verification scripts by @pauldruce in [#75](https://github.com/pauldruce/RFL/pull/75)

### 📚 Documentation & Governance
* docs(roadmap): establish v0.3.0 release roadmap, issue tracking, and EP alignment by @pauldruce in [#72](https://github.com/pauldruce/RFL/pull/72)
* docs(eps): propose EP-4 computational physics verification suite by @pauldruce in [#64](https://github.com/pauldruce/RFL/pull/64)

## [v0.2.0] - 2026-09-13

### Highlights & Breaking Changes
* **Canonical Root Packaging & Downstream Isolation:** Relocates build manifests and directories to the repository root. Isolates internal test targets and GoogleTest from downstream `FetchContent` consumers ([#54](https://github.com/pauldruce/RFL/issues/54), [#55](https://github.com/pauldruce/RFL/pull/55)).
* **AddressSanitizer Gating & Optimisation:** Replaces unconditional ASan flags with `RFL_ENABLE_ASAN`, eliminating symbol pollution in static libraries and accelerating test execution ([#47](https://github.com/pauldruce/RFL/issues/47), [#55](https://github.com/pauldruce/RFL/pull/55)).
* **Streamlined Build Architecture & Please Removal:** Removes obsolete Please build manifests in favour of CMake, `scikit-build-core`, and `mise` ([#55](https://github.com/pauldruce/RFL/pull/55)).
* **Expanded Example Workflows:** Adds CMake integration, standalone `Makefile` targets, and direct `g++` compilation workflows across example directories ([#55](https://github.com/pauldruce/RFL/pull/55)).
* **Native Windows MSVC Support:** Compiles natively under MSVC on Windows and publishes precompiled `win_amd64` wheels to PyPI ([#21](https://github.com/pauldruce/RFL/issues/21), [#48](https://github.com/pauldruce/RFL/pull/48), [#50](https://github.com/pauldruce/RFL/pull/50)).
* **Dependency Decoupling with `FetchContent`:** CMake automatically fetches Armadillo when missing from the system ([#22](https://github.com/pauldruce/RFL/issues/22), [#43](https://github.com/pauldruce/RFL/pull/43)).
* **Resource-Efficient CI Pipeline:** Introduces job-level path filtering, compiler caching (`ccache`), and tiered PR smoke tests to cut compute time ([#42](https://github.com/pauldruce/RFL/pull/42)).
* **Release Automation & Integrity Gates:** Enforces dynamic Git tag versioning, PR milestone verification, documentation link checking, and changelog release gates ([#51](https://github.com/pauldruce/RFL/pull/51)).
* **Breaking Change (Python 3.8):** Dropped Python 3.8 support to align with `cibuildwheel` v4 and active scientific Python baselines ([#40](https://github.com/pauldruce/RFL/issues/40)).

### 🚀 Features & Enhancements
* feat(examples): add CMake and standalone Makefile workflows with comprehensive documentation by @pauldruce in [#55](https://github.com/pauldruce/RFL/pull/55)
* feat(ci,cmake): implement native Windows support and Python wheels ([#21](https://github.com/pauldruce/RFL/issues/21)) by @pauldruce in [#48](https://github.com/pauldruce/RFL/pull/48)
* feat(cmake,ci): decouple dependencies using FetchContent and unify CI workflows ([#22](https://github.com/pauldruce/RFL/issues/22)) by @pauldruce in [#43](https://github.com/pauldruce/RFL/pull/43)
* feat(ci,release): dynamic versioning, changelog automation, and release integrity gates by @pauldruce in [#51](https://github.com/pauldruce/RFL/pull/51)
* feat(packaging): dynamic PEP 440 pre-release versioning for PyPI release candidates by @pauldruce in [#41](https://github.com/pauldruce/RFL/pull/41)

### 🐛 Bug Fixes
* fix(cmake): migrate CMakeLists.txt and pyproject.toml to root, standardise source layout in src/core, and isolate downstream consumer test targets ([#54](https://github.com/pauldruce/RFL/issues/54)) by @pauldruce in [#55](https://github.com/pauldruce/RFL/pull/55)
* fix(cmake): gate AddressSanitizer behind RFL_ENABLE_ASAN to prevent static library symbol leakage ([#47](https://github.com/pauldruce/RFL/issues/47)) by @pauldruce in [#55](https://github.com/pauldruce/RFL/pull/55)
* fix(ci,windows): install lapack in vcpkg and isolate armadillo cache keys by @pauldruce in [#50](https://github.com/pauldruce/RFL/pull/50)
* fix(ci,cmake): support legacy dependencies under CMake 4.x and declutter step summary by @pauldruce in [#44](https://github.com/pauldruce/RFL/pull/44)
* fix(ci): bump actions to Node 24 runtimes to eliminate deprecation warnings by @pauldruce in [#45](https://github.com/pauldruce/RFL/pull/45)

### 🧰 Build & CI/CD Architecture
* release(tooling): establish pre-release qualification suite, isolate example outputs, and harden changelog verification by @pauldruce in [#63](https://github.com/pauldruce/RFL/pull/63)
* chore(build): remove obsolete Please build system and enforce repository-wide CSpell checks by @pauldruce in [#55](https://github.com/pauldruce/RFL/pull/55)
* fix(ci,release): improve changelog PR detection and self-reference release PRs by @pauldruce in [#53](https://github.com/pauldruce/RFL/pull/53)
* chore(release): prepare v0.2.0 release and CHANGELOG.md by @pauldruce in [#52](https://github.com/pauldruce/RFL/pull/52)
* ci(deps): bump actions/setup-python from 5 to 7 in the github-actions group across 1 directory by @dependabot[bot] in [#61](https://github.com/pauldruce/RFL/pull/61)
* ci(deps): bump the github-actions group with 7 updates by @dependabot[bot] in [#37](https://github.com/pauldruce/RFL/pull/37)

### 📚 Documentation & Governance
* docs: document example workflows, migrate active backlog tasks to GitHub issues, and remove obsolete v0.1.0 release document ([#57](https://github.com/pauldruce/RFL/issues/57), [#58](https://github.com/pauldruce/RFL/issues/58), [#59](https://github.com/pauldruce/RFL/issues/59), [#60](https://github.com/pauldruce/RFL/issues/60)) by @pauldruce in [#55](https://github.com/pauldruce/RFL/pull/55)
* docs(eps): add EP-3 for resource-efficient CI/CD pipeline & modern build architecture by @pauldruce in [#42](https://github.com/pauldruce/RFL/pull/42)
* docs: post-v0.1.0 release updates, README overhaul, and PyPI metadata by @pauldruce in [#39](https://github.com/pauldruce/RFL/pull/39)

**Full Changelog**: https://github.com/pauldruce/RFL/compare/v0.1.0...v0.2.0

## [v0.1.0] - 2026-08-31

### Highlights
* Initial packaged release of Random Fuzzy Library (`pyrfl` on PyPI).
* Automated multi-platform Python wheel packaging for Linux (`manylinux_2_28`) and macOS (`x86_64`, `arm64`) across Python 3.9–3.13.
* Zero-copy NumPy array interoperability using Carma and pybind11.
* Restructured codebase into `src/RFL/core` (`RFL::core`) and `src/RFL/legacy` (`RFL::legacy`).
* Established the Enhancement Proposal (EP) framework with EP-1 and EP-2.
* Standardised technical documentation and C++ docstrings using ASD-STE100 and British English.

[v0.3.0]: https://github.com/pauldruce/RFL/compare/v0.2.0...HEAD
[v0.2.0]: https://github.com/pauldruce/RFL/compare/v0.1.0...v0.2.0
[v0.1.0]: https://github.com/pauldruce/RFL/releases/tag/v0.1.0
