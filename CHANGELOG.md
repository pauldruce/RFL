# Changelog

All notable changes to the Random Fuzzy Library (RFL) are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
* chore(build): remove obsolete Please build system and enforce repository-wide CSpell checks by @pauldruce in [#55](https://github.com/pauldruce/RFL/pull/55)
* fix(ci,release): improve changelog PR detection and self-reference release PRs by @pauldruce in [#53](https://github.com/pauldruce/RFL/pull/53)
* chore(release): prepare v0.2.0 release and CHANGELOG.md by @pauldruce in [#52](https://github.com/pauldruce/RFL/pull/52)
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
