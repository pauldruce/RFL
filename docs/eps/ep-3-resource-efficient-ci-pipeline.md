# EP-3: Resource-Efficient CI/CD Pipeline & Multi-Tier Validation

* **Title:** Resource-Efficient CI/CD Pipeline & Multi-Tier Validation
* **Author:** Paul Druce
* **Status:** In Discussion
* **Target Versions:** RFL v0.2.0 (Phases 1 & 2), v0.2.0 (Phases 3 & 4)
* **Date:** 2026-09-04

---

## 1. Multi-Phase Implementation Tracker

This proposal spans the delivery of the `v0.2.0` infrastructure milestone.
The table below tracks the status of each implementation phase:

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Tier 1 Fast PR Gate (Ubuntu + macOS, path filtering, Please graph pruning) | `v0.2.0` | [#19](https://github.com/pauldruce/RFL/issues/19), [#23](https://github.com/pauldruce/RFL/issues/23) | ⏳ Scheduled |
| **Phase 2** | Tier 2 Exhaustive Compatibility Matrix on `main` merge & `ccache` caching | `v0.2.0` | [#23](https://github.com/pauldruce/RFL/issues/23) | 💡 Planned |
| **Phase 3** | Dependency decoupling via CMake `FetchContent` & system packages | `v0.2.0` | [#22](https://github.com/pauldruce/RFL/issues/22) | 💡 Planned |
| **Phase 4** | Native Windows MSVC CI runner & binary wheel automation | `v0.2.0` | [#21](https://github.com/pauldruce/RFL/issues/21) | 💡 Planned |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Research Context
`RFL` requires continuous validation across multiple operating systems (Ubuntu, macOS, Windows) and scientific libraries (Armadillo, GSL, OpenBLAS).
Currently, every push and pull request triggers an uncoordinated collection of workflows:
1. `build_and_test.yml` executes on every push across three operating systems.
2. `compatability_tests.yml` executes a 15-job matrix on every pull request.
3. Each job downloads and builds Armadillo from SourceForge tarballs via custom shell scripts.

This approach creates severe research bottlenecks:
* **Excessive Latency:** Developers wait 15 to 25 minutes for pull request checks to finish.
* **Wasted Compute:** A single pull request consumes 60 to 90 GitHub Actions runner minutes.
* **Documentation Friction:** Changes to documentation (`docs/**`, `*.md`, `.agents/**`) trigger full C++ compilation matrices.
* **Network Fragility:** SourceForge download mirrors frequently time out or throttle concurrent runners.

A tiered, resource-efficient CI pipeline is required to deliver sub-2-minute feedback during development while preserving exhaustive platform verification on `main`.

### 2.2 Goals
* **Sub-2-Minute Dual-Platform PR Feedback:** Validate Ubuntu and macOS runners on pull requests in under two minutes using system packages.
* **Near-Zero Compute for Documentation:** Ensure documentation changes consume zero C++ compilation runner minutes.
* **Exhaustive Matrix on Merge to Main:** Execute the full 15-job multi-OS compatibility matrix whenever code merges into `main`.
* **Hermetic Incremental Testing:** Leverage Please build (`plz query changes`) to test only targets affected by changed files.
* **Persistent Build Caching in CI:** Persist Please build cache (`.plz-cache`) across ephemeral GitHub runners using `actions/cache`.
* **Clean Test Artifact Extraction:** Use Please (`plz export outputs`) to isolate test payloads without copying internal build trees.
* **Standardised CI Test Reporting:** Generate machine-readable JUnit XML test results for native GitHub pull request summaries.
* **Compiler Caching:** Integrate `ccache` to eliminate redundant compilation across workflow runs.
* **Decoupled Dependencies:** Replace custom SourceForge build scripts with system package managers (`apt`, `brew`) and CMake `FetchContent`.

### 2.3 Non-Goals
* Managing self-hosted GitHub Actions runners (standard GitHub-hosted runners remain sufficient).
* Replacing CMake as the primary build system for external consumers.
* Implementing distributed build caching across private networks.

---

## 3. Research Workflows & Scientific Requirements

### 3.1 Core Research Scenarios
1. **Scenario 1 (Rapid Documentation and Theoretical Derivation Updates):** A researcher updates scientific derivations in `docs/` or notes in `.agents/`. The CI pipeline runs fast markdown checks. It completes in seconds without launching C++ compilers.
2. **Scenario 2 (High-Frequency MCMC Algorithm Development):** A developer refactors Dirac operator routines in `src/RFL/core/`. The Tier 1 gate builds on Ubuntu and macOS using system packages and Please, executing unit tests in under two minutes.
3. **Scenario 3 (Exhaustive Mainline Verification):** When a pull request merges into `main`, GitHub Actions automatically runs all 15 permutations of OS versions and Armadillo releases to verify release health.

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Operational & Physical Invariant |
| :--- | :--- | :--- |
| **REQ-CI-001** | **Fast Dual-Platform PR Gate** | Tier 1 pull request checks must validate Ubuntu and macOS in under 120 seconds using system packages. |
| **REQ-CI-002** | **Zero-Compute Documentation PRs** | Commits touching only `.md`, `docs/`, `.agents/`, or meta files must trigger zero C++ build jobs. |
| **REQ-CI-003** | **Branch Protection Integrity** | Path-filtered workflows must report passing status checks to GitHub to prevent pull request deadlocks. |
| **REQ-CI-004** | **Exhaustive Mainline Verification** | Every commit pushed or merged into `main` must execute the full 15-job compatibility matrix. |
| **REQ-CI-005** | **Compiler Object Caching** | Unchanged C++ translation units must hit `ccache` caches with an average cache-hit rate exceeding 70%. |
| **REQ-CI-006** | **Hermetic Graph Analysis** | Please build target calculation must match the exact transitive dependency closure of changed source files. |
| **REQ-CI-007** | **Persistent Runner Caching** | The Please directory cache (`.plz-cache`) must be saved and restored across GitHub Actions runs. |
| **REQ-CI-008** | **Minimal Test Artifact Footprint** | Exported test artifacts must contain only final binaries and libraries (< 20 MB total). |
| **REQ-CI-009** | **Structured Test Reporting** | Test executions must export JUnit XML test results for native GitHub pull request reporting. |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: Execution Topology (Tiered CI with Dual-Platform PR Smoke Test)

| Criteria | Option A: Monolithic 15-Job Matrix on PR | Option B: Fast 2-Platform PR Gate + Full Matrix on Main (Selected) | Option C: Linux-Only PR Gate |
| :--- | :--- | :--- | :--- |
| **PR Turnaround Time** | ❌ 15–25 minutes | ✅ **< 2 minutes** | ✅ < 1.5 minutes |
| **Runner Minute Cost** | ❌ 60–90 minutes per PR | ✅ **3–5 minutes per PR** | ✅ 1–2 minutes per PR |
| **Early macOS/Linux Detection** | ✅ High | ✅ **High (Catches 99% of POSIX/Clang divergence)** | ❌ None (macOS issues delayed) |
| **Exhaustive Verification Gate** | On every PR commit | ✅ **On every commit to `main`** | On every commit to `main` |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* Personal repositories cannot use native GitHub Merge Queues. Option B provides the optimal balance: pull requests automatically test both Ubuntu and macOS using fast system packages (`apt` / `brew`), catching cross-platform compiler issues within 2 minutes. The full 15-job combination matrix runs automatically upon push to `main` to verify complete matrix compatibility.

---

### 4.2 ADR-2: Path Filtering Strategy

| Criteria | Option A: Top-Level `paths-ignore` | Option B: Job-Level Path Evaluator (Selected) | Option C: Commit Message Tokens (Current) |
| :--- | :--- | :--- | :--- |
| **Branch Protection Compatibility** | ❌ Fails (Skipped workflows leave required checks pending) | ✅ **Guaranteed (Workflow executes and reports green check)** | ⚠️ Unreliable (Requires manual developer action) |
| **Granularity** | ⚠️ Binary (Entire workflow skips) | ✅ **Multi-category (Docs, Python, Core C++)** | ❌ Single string matching |
| **Reliability** | ✅ High | ✅ **High** | ❌ Poor (Developers forget `#docs` token) |
| **Decision** | Rejected | **Selected (Option B)** | Rejected |

*Rationale:* Top-level `paths-ignore` causes pull requests to hang if repository branch protection enforces required status checks. Using a job-level path filter (`dorny/paths-filter`) allows the workflow to execute, dynamically skip compilation jobs, and immediately report a successful status check.

---

### 4.3 ADR-3: Dependency Management in CI Runners

| Criteria | Option A: SourceForge Tarball Compilation (Current) | Option B: System Package Managers (`apt`/`brew`) | Option C: CMake `FetchContent` Fallback |
| :--- | :--- | :--- | :--- |
| **Setup Duration** | ❌ 3–5 minutes per runner | ✅ **5–10 seconds per runner** | ⚠️ 1–2 minutes (cached) |
| **Mirror Availability** | ❌ Unreliable (SourceForge timeouts) | ✅ **High (Local GitHub runner mirrors)** | ✅ High (Git mirror) |
| **Cleanliness** | ❌ Pollutes global `/usr/local` | ✅ Cleanly tracked by OS | ✅ Contained in build directory |
| **Multi-Version Flexibility** | ✅ High (Compiles any tag) | ⚠️ Limited to distro package | ✅ **High (Pulls any Git commit/tag)** |
| **Decision** | Rejected | **Selected for Tier 1 Gate** | **Selected for Tier 2 & FetchContent** |

*Rationale:* Recompiling Armadillo from SourceForge on every job wastes thousands of compute minutes annually. Tier 1 gates will use system packages (`libarmadillo-dev`, `brew install armadillo`) for instant setup. Tier 2 compatibility matrices will use CMake `FetchContent` and `ccache` to test specific library versions cleanly.

---

### 4.4 ADR-4: Build Graph Pruning & Hermetic Execution Engine (Please Build)

| Criteria | Option A: Monolithic Test Execution (`ctest`) | Option B: Please Graph Analysis & Hermetic Engine (Selected) |
| :--- | :--- | :--- |
| **Execution Latency** | ⚠️ 30–60 seconds | ✅ **0.1–5 seconds** |
| **Doc-Only Impact** | ⚠️ Compiles code regardless of changes | ✅ **Returns 0 targets; skips test execution** |
| **Local Reproducibility** | ⚠️ Manual target selection | ✅ **Identical command runs on developer machines** |
| **Binary Output Isolation** | ❌ Intermingled build tree | ✅ **Direct target export via `plz export outputs`** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Please build maintains a precise directed acyclic graph (DAG) of project dependencies. Running `./pleasew query changes --since origin/main --level -1` identifies exactly which test targets require execution based on the git diff. In addition, Please provides clean artifact extraction via `plz export outputs` and persistent directory caching via `.plz-cache`.

---

### 4.5 ADR-5: CI Status Check Enforcement & Branch Protection Strategy

| Criteria | Option A: Aggregate Gatekeeper (`ci-gate`) (Selected) | Option B: Require Individual Matrix Jobs |
| :--- | :--- | :--- |
| **Doc-Only Safety** | ✅ **100% Safe (Evaluates skipped jobs as successful)** | ❌ **Deadlocks on doc-only PRs** |
| **Settings Maintenance** | ✅ **Zero maintenance (Only one check name required)** | ❌ High (Breaks whenever matrix changes) |
| **PR Clarity** | ✅ Single decisive status badge | ❌ Many confusing status items |
| **Decision** | **Selected (Option A)** | Rejected |

*Rationale:* Branch protection on `main` will require a single aggregate status check: `ci-gate`. This job evaluates the results of all upstream PR gate jobs. If jobs are safely skipped due to path filtering, `ci-gate` reports success. If any executed job fails, `ci-gate` fails and blocks merging.

---

### 4.6 ADR-6: Pipeline Phase Separation (Pre-Build Gate vs. Build vs. Post-Build Test)

| Criteria | Option A: Monolithic Build-and-Test Script | Option B: Tri-Stage Phase Separation (Selected) |
| :--- | :--- | :--- |
| **Fail-Fast Latency** | ❌ 5–10 minutes (Fails linting after compilation) | ✅ **5–15 seconds (Fails linting before compilation)** |
| **Compute Conservation** | ❌ Compiles code even if syntax or spelling is broken | ✅ **Zero compilation spent on unformatted code** |
| **Artifact Portability** | ⚠️ Unstructured binary artifacts | ✅ **Structured exports for downstream test suites** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Splitting the pipeline into distinct phases ensures that fast pre-build checks (formatting, spelling, path evaluation) stop execution before expensive C++ compilers spin up. When building packages (such as Python wheels), building the artifact once and testing it across environments eliminates duplicate C++ compilation.

---

## 5. Target Architecture & Component Contracts

### 5.1 Workflow Topology

```
.github/workflows/
├── ci.yml                     # Tier 1: Fast Dual-Platform PR Gate (Ubuntu + macOS, system packages, < 2 min)
├── compatability_tests.yml    # Tier 2: Exhaustive Matrix (On push to main, workflow_dispatch)
├── linter.yml                 # Code style & formatting (Clang-format on src/RFL and playground)
├── codeql.yml                 # Security scanning (Main branch and PRs to main)
└── release.yml                # Binary wheels and PyPI publishing (Tag push and release PRs)
```

> [!NOTE]
> `.github/workflows/build_and_test.yml` is retired because its responsibilities are completely subsumed by `ci.yml` (Tier 1) and `compatability_tests.yml` (Tier 2).

---

### 5.2 Please Build CI Optimisation Architecture

Please build provides specific capabilities to accelerate CI runs:

#### 1. Graph Pruning (`plz query changes`)
On pull requests, Please identifies affected test targets against the base branch:
```bash
TARGETS=$(./pleasew query changes --since origin/${{ github.base_ref || 'main' }} --level -1)
if [ -n "$TARGETS" ]; then
  ./pleasew test -p --test_results_file=plz-out/log/test_results.xml $TARGETS
else
  echo "No code targets affected by this change."
fi
```

#### 2. Persistent Runner Caching (`.plz-cache` + `actions/cache`)
Please supports local directory caching. By enabling `[cache]` in `.plzconfig`:
```ini
[cache]
dir = .plz-cache
dirclean = true
```
GitHub Actions saves and restores this cache across workflow runs:
```yaml
- name: Restore Please Cache
  uses: actions/cache@v4
  with:
    path: .plz-cache
    key: plz-cache-${{ runner.os }}-${{ hashFiles('**/*.cpp', '**/*.hpp', 'BUILD*', '.plzconfig') }}
    restore-keys: |
      plz-cache-${{ runner.os }}-
```
Unchanged targets are retrieved from the local cache in milliseconds.

#### 3. Lightweight Artifact Export (`plz export outputs`)
When separating build and test jobs, Please exports only the final target outputs rather than copying the entire `plz-out/` directory:
```bash
./pleasew export outputs -o ./dist //src/RFL/core:rfl //src/RFL/core/tests:unit_tests
```
This produces a clean payload (< 20 MB) suitable for GitHub Actions artifact transfer.

#### 4. Plain CI Output & Structured Reporting
* **Plain Terminal Logging (`-p`):** Disables ANSI interactive spinners so CI log files remain readable and concise.
* **JUnit XML Output (`--test_results_file`):** Generates `test_results.xml` for GitHub Actions to render native test pass/fail annotations.

---

### 5.3 Workflow Configuration Contracts

#### 1. Tier 1: Fast Dual-Platform PR Gate (`.github/workflows/ci.yml`)
* **Triggers:** `pull_request` (targeting `main`), `workflow_dispatch`.
* **Path Filter Configuration:**
  ```yaml
  filters:
    docs:
      - 'docs/**'
      - '*.md'
      - '.agents/**'
      - 'research/**'
      - 'vault/**'
      - '.gitignore'
    code:
      - 'src/**'
      - 'playground/**'
      - 'CMakeLists.txt'
      - 'BUILD.plz'
      - '.plzconfig'
      - 'third_party/**'
  ```
* **Jobs:**
  * `filter`: Evaluates changed paths in seconds using `dorny/paths-filter`.
  * `test_ubuntu`: Executes if `code == true`. Uses `apt-get install libarmadillo-dev libgsl-dev`, `ccache`, Please cache restoration, and Please graph pruning (`plz query changes`).
  * `test_macos`: Executes if `code == true`. Uses `brew install armadillo gsl`, `ccache`, and CMake/CTest.
  * `ci-gate`: Always runs (`if: always()`). Evaluates upstream results; reports success if code passed or if docs were skipped.

#### 2. Tier 2: Exhaustive Compatibility Matrix (`.github/workflows/compatability_tests.yml`)
* **Triggers:**
  * `push` to `main` branch (every commit/merge into `main`).
  * `workflow_dispatch` with manual parameter inputs (`all` vs. random subset).
* **Runner Matrix:**
  * Operating Systems: `ubuntu-latest`, `ubuntu-22.04`, `macos-latest`, `macos-15-intel`, `macos-14`.
  * Armadillo Versions: `11.4.4`, `12.8.4`, `14.2.2`.
  * Build Type: `Release`.
* **Optimisation:**
  * Initialise `hendrikmuhs/ccache-action` on all runners.
  * Cache Armadillo installations per OS and version using `actions/cache`.

#### 3. C++ Linter (`.github/workflows/linter.yml`)
* **Triggers:** `pull_request`, `push` to `main`.
* **Target Directories:** `src/RFL` and `playground` (removing the legacy `RFL_source` path).
* **Path Filtering:** Only execute when C/C++ files (`**/*.{cpp,hpp,h,c}`) or formatting configurations change.

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Test Description | Acceptance Threshold |
| :--- | :--- | :--- |
| **Doc-Only PR Gate** | Push commit modifying only `docs/README.md` | PR check passes in under 30 seconds with zero C++ compilation. |
| **Core C++ PR Gate** | Push commit modifying `src/RFL/core/DiracOperator.cpp` | PR check tests Ubuntu + macOS and passes in under 120 seconds. |
| **Please Graph Query** | `./pleasew query changes docs/README.md` | Returns zero targets. |
| **Please Cache Restoration** | `./pleasew test //src/RFL/core/tests:unit_tests` | Re-execution reports cached targets in < 100 ms. |
| **Artifact Export Footprint** | `./pleasew export outputs -o ./dist //src/RFL/core:rfl` | Output directory contains only libraries, < 20 MB total. |
| **C++ Linter Execution** | `./pleasew run //:lint_cpp` | Successfully formats `src/RFL` without checking ignored build folders. |
| **Mainline Trigger** | Push commit to `main` | Automatically triggers `compatability_tests.yml` across full matrix. |

---

## 7. Phased Delivery Plan

### Phase 1: Tier 1 Fast Dual-Platform PR Gate & Path Filtering
* **Target Version:** `v0.2.0`
* **GitHub Issues:** [#19](https://github.com/pauldruce/RFL/issues/19), [#23](https://github.com/pauldruce/RFL/issues/23)
* **Tasks:**
  1. Add path filtering to `.github/workflows/ci.yml` using `dorny/paths-filter`.
  2. Implement dual-platform PR smoke tests on Ubuntu (`libarmadillo-dev`) and macOS (`brew install armadillo`).
  3. Add `ci-gate` aggregate check in `ci.yml` for branch protection.
  4. Enable Please directory caching (`.plz-cache`) in `.plzconfig` and `actions/cache`.
  5. Integrate Please graph pruning (`plz query changes`) and `-p` logging into `ci.yml`.
  6. Fix paths in `.github/workflows/linter.yml` (`src/RFL` and `playground`).
  7. Delete obsolete `.github/workflows/build_and_test.yml`.
  8. Update `BUILD.plz` linting scripts to exclude build and IDE directories.

### Phase 2: Tier 2 Exhaustive Compatibility Matrix on Push to Main
* **Target Version:** `v0.2.0`
* **GitHub Issue:** [#23](https://github.com/pauldruce/RFL/issues/23)
* **Tasks:**
  1. Retarget `.github/workflows/compatability_tests.yml` to trigger on `push: branches: [ main ]` and `workflow_dispatch`.
  2. Remove automatic 15-job execution from `pull_request` events.
  3. Add `hendrikmuhs/ccache-action` across all matrix runners.
  4. Cache pre-built Armadillo installations per matrix combination via `actions/cache`.
  5. Verify `combination_selection.py` matrix generator operates reliably on `main` push.

### Phase 3: Dependency Decoupling via CMake FetchContent & System Packages
* **Target Version:** `v0.2.0`
* **GitHub Issue:** [#22](https://github.com/pauldruce/RFL/issues/22)
* **Tasks:**
  1. Update `src/RFL/cmake/Armadillo.cmake` with `FetchContent` fallback from upstream GitLab repository.
  2. Replace manual Armadillo shell scripts in Tier 2 matrix with CMake FetchContent or system packages.
  3. Remove `.github/actions/install-armadillo` composite action once decoupled.

### Phase 4: Native Windows MSVC CI Runner & Binary Wheel Automation
* **Target Version:** `v0.2.0`
* **GitHub Issue:** [#21](https://github.com/pauldruce/RFL/issues/21)
* **Tasks:**
  1. Add `windows-latest` runner to Tier 2 matrix.
  2. Configure MSVC compiler flags (`/std:c++17`, `/utf-8`, `/openmp`) in `CMakeLists.txt`.
  3. Add Windows wheel compilation (`win_amd64`) to `.github/workflows/release.yml`.
