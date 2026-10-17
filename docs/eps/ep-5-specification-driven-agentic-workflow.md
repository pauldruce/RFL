# EP-5: Specification-Driven Agentic Workflow and Traceability Infrastructure

* **Title:** Specification-Driven Agentic Workflow and Traceability Infrastructure
* **Author:** Paul Druce
* **Status:** Accepted
* **Target Version:** RFL v0.3.0
* **Date:** 2026-10-14

---

## 1. Multi-Phase Implementation Tracker

This proposal establishes a lightweight workflow and traceability infrastructure for human architects and autonomous AI coding agents.
The table below tracks the status of each implementation phase:

| Phase | Scope and Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Global specification schema and requirement namespacing standard (`REQ-<SYS>-<NUM>`) across EPs | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Implemented |
| **Phase 2** | Zero-dependency requirement-to-test traceability engine across C++ and Python test suites (ADR-1) | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Implemented |
| **Phase 3** | Automated change impact analysis and selective test execution planner across C++ and Python (ADR-3) | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Implemented |
| **Phase 4** | Shift-left developer feedback loops and pull request quality gate summaries via encapsulated task contracts (ADR-2, ADR-4) | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Implemented |

---

## 2. Motivation, Goals and Non-Goals

### 2.1 Problem Statement and Agentic Context
Software development in `RFL` pairs a human physicist-architect with autonomous AI coding agents.
Without clear structural guardrails, agentic workflows exhibit four recurring failure modes:

1. **Requirement Collision and Ambiguity:** Un-namespaced identifiers like `REQ-001` in multiple EPs prevent deterministic lookup across specifications and test suites.
2. **Traceability Drift:** Scientific invariants defined in EPs drift from test suites because no automated tool links specifications to executable tests.
3. **Slow Feedback Loops:** Relying on remote CI pipelines adds a multi-minute delay per iteration. Agents finish active turns before discovering defects.
4. **Coarse Invalidation:** Static file path matching executes redundant test suites or misses indirect cross-language dependencies between C++ and Python.

### 2.2 Goals
* **Deterministic Verification:** Provide zero-dependency automated verification that maps requirement specifications (`REQ-<SYS>-<NUM>`) to test cases across C++ and Python test frameworks.
* **Domain-Aware Auditing:** Distinguish unit-testable code and physics invariants (`ARCH`, `PHY`, `TRC`) from CI and packaging rules (`CI`, `PKG`).
* **Selective Execution Planning:** Deterministically identify affected tests and impacted requirements from code and specification changes without requiring full regression runs.
* **Shift-Left Agentic Feedback:** Execute fast static checks during local pre-commit in under 1 second, and run affected test suites before push.
* **Low-Friction PR Governance:** Publish automated traceability tables and modified-test diff alerts to pull request summaries without blocking refactoring PRs.
* **Preserve CI Budget:** Complete all traceability audits and execution planning in under 200 ms, preserving the sub-120-second CI budget (`REQ-CI-001`).

### 2.3 Non-Goals
* Do not adopt heavyweight external ALM frameworks such as Sphinx-Needs, OpenFastTrace, ReqIF, or Doorstop.
* Do not enforce rigid multi-commit Git history choreography that breaks during interactive rebasing or debugging.
* Do not require cryptographic SHA-256 manifest files or special GitHub override labels to modify test suites.
* Do not track auto-generated traceability matrix files in Git where they cause merge conflicts.
* Do not pollute production source code in `src/` with requirement annotations or CI directives.

---

## 3. Requirements and Operational Invariants

| Requirement ID | Requirement Summary | Lifecycle Status | Operational and Architectural Invariant |
| :--- | :--- | :--- | :--- |
| **REQ-TRC-001** | **Global Requirement Namespacing** | Phase 1 (Implemented) | Requirement IDs must conform to regex `^REQ-[A-Z]{2,4}-\d{3}$`. Subsystem codes: `ARCH` (Architecture), `PKG` (Packaging), `CI` (CI/CD), `PHY` (Physics), `TRC` (Traceability). IDs must remain globally unique across all EPs. |
| **REQ-TRC-002** | **Lifecycle and Domain Awareness** | Phase 1 (Implemented) | Section 3 requirement tables across all EPs must maintain four columns: `Requirement ID`, `Requirement Summary`, `Lifecycle Status`, and an invariant definition column (`Operational and Architectural Invariant`). Traceability audits evaluate unit-testable requirements (`ARCH`, `PHY`, `TRC`) marked `Implemented`. |
| **REQ-TRC-003** | **Zero-Dependency Traceability Engine** | Phase 2 (Implemented) | The traceability verification engine must use Python ($\ge 3.10$) standard library modules only. Process exit codes must be deterministic: 0 (pass), 1 (validation error), 2 (command-line usage error). Execution must complete in under 200 ms. |
| **REQ-TRC-004** | **Multi-Framework Tag Support** | Phase 2 (Implemented) | The verification engine must extract tags from GoogleTest comments (`// [REQ-...]`), Catch2 tags (`[REQ-...]`), and Pytest test functions (`def test_*` via marks, docstrings, or comments), ignoring non-test helper methods. |
| **REQ-TRC-005** | **Automated Change Impact Analysis** | Phase 3 (Implemented) | Given a changeset across source files, bindings, build configurations, or specifications, the execution planner must deterministically identify all affected test suites and impacted requirements without source code annotations. Verification gates must validate test results via JUnit XML and highlight modified test files in pull request summaries. |
| **REQ-TRC-006** | **Shift-Left Local Git Hooks** | Phase 4 (Implemented) | Local git hooks must enforce fast pre-commit static validation (spelling, documentation links, specification audits) in under 1 second. Pre-push validation must run linting and affected component and requirement test suites to block broken changes before remote transmission. |

---

## 4. Architecture Decision Records (ADRs) and Trade-offs

### 4.1 ADR-1: Zero-Dependency Python Tooling vs External ALM Frameworks

| Criteria | Option A: External Tools (Sphinx-Needs / Doorstop) | Option B: Custom Standard Library Python Script (Selected) |
| :--- | :--- | :--- |
| **Dependencies** | Heavy (Requires Sphinx build or external Python packages) | **Zero (Standard Python 3 only, native to RFL)** |
| **Execution Latency** | 5–30 seconds in CI | **< 200 ms (instant local and CI feedback)** |
| **Specification Format** | Poor (Requires YAML files or rST directives) | **Native (Parses GitHub Flavoured Markdown tables directly)** |
| **Maintenance** | Vulnerable to upstream deprecation and package churn | **Complete repository ownership co-located in `scripts/`** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* RFL prioritises lean engineering principles. Heavy ALM frameworks violate `REQ-CI-001` and `REQ-CI-005` by inflating CI runtimes and adding external dependencies. A custom Python 3 script (`scripts/verify_traceability.py`) uses the standard library only. It parses GitHub Flavoured Markdown tables directly, executes in under 200 ms, and integrates cleanly with `mise`.

---

### 4.2 ADR-2: JUnit Proof and Test-Diff Visibility vs Git Choreography and SHA-256 Manifest Ceremonies

| Criteria | Option A: Two-Commit Git Red-Green Choreography & SHA-256 Manifests | Option B: JUnit XML Verification & PR Test-Diff Visibility (Selected) |
| :--- | :--- | :--- |
| **C++ Iteration Ergonomics** | ❌ Poor (Fights interactive debugging, rebasing, and fixup commits) | ✅ **Optimal (Standard atomic commits across headers, source, and tests)** |
| **Test Addition Friction** | ❌ High (Editing tests or fixing tolerances breaks hash manifests) | ✅ **Zero (New or modified tests require no manifest regeneration)** |
| **Tamper Visibility** | Cryptographic hash mismatch requires manual override label | **Explicit alert in** `$GITHUB_STEP_SUMMARY` **listing modified test files** |
| **Execution Proof** | Fragile Git checkout time-travel in CI | **Direct verification of CTest and Pytest JUnit XML results** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Enforcing multi-commit Git sequences or SHA-256 test manifests introduces high rigidity without preventing subtle test tampering. Autonomous agents often create trivial passing tests to satisfy commit-order checks. In RFL, a human architect reviews every pull request. Combining automated JUnit XML test proof with an explicit test-diff alert in `$GITHUB_STEP_SUMMARY` provides robust verification without ceremonial overhead.

---

### 4.3 ADR-3: Automated Dependency Impact Planning vs Source Code Tagging and Static Path Filters

| Criteria | Option A: Source Code Tagging & Static `dorny/paths-filter` | Option B: Automated Dependency Graph via `#include` and `pybind11` AST (Selected) |
| :--- | :--- | :--- |
| **Source Code Cleanliness** | ❌ Poor (Pollutes C++ headers and sources with governance tags) | ✅ **100% Clean (`src/` contains zero requirement tags or CI comments)** |
| **Cross-Language Semantic Bridge** | Manual `# +ci:depends_on` comments that rot over time | **Automatic (`bindings.cpp` symbol table + Python AST `rfl.<Symbol>`)** |
| **Specification Edit Handling** | Static path filters skip tests when only docs change | **Automatic (`git diff` on EP tables schedules mapped C++ and Python tests)** |
| **CI Execution Speed** | Coarse file filters trigger redundant full-suite test runs | **Runs only affected CTest suites (`-R`) and affected Python test files** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Production code in `src/` must remain pure C++ and Python without requirement annotations or magic CI comments. Static path filters like `dorny/paths-filter` lack semantic code awareness and cannot trace cross-language dependencies. Unit tests sit between the production code they exercise and the `[REQ-...]` tags they verify. The planner (`scripts/ci_plan.py`) resolves impact across two hops:
1. **Hop 1 (`Code -> Affected Tests`):** Traces C++ `#include` graphs and resolves Python test ASTs (`rfl.<Symbol>`) through `bindings.cpp` (`py::class_`).
2. **Hop 2 (`Affected Tests + EP Diffs -> Impacted Requirements`):** Extracts `[REQ-...]` tags from affected tests, modified EP table rows, and `Implements: REQ-...` Git trailers.
3. **Hop 3 (`Impacted Requirements -> Full Requirement Test Closure`):** Expands impacted requirements to their full test closure, emitting `plan.json` for targeted CTest (`-R`) and Pytest execution.

---

### 4.4 ADR-4: Shift-Left Local Git Hooks vs Remote-Only CI Gating

| Criteria | Option A: Remote CI Only | Option B: Shift-Left Dual Loop (Local Hooks + Remote CI) (Selected) |
| :--- | :--- | :--- |
| **Agent Feedback Latency** | 2–5 minutes per push | **< 1 second (instant stderr self-correction in local turn)** |
| **PR Commit Hygiene** | Poor (noisy retry commits in pull request history) | **Optimal (clean commits verified before push)** |
| **Compute Consumption** | High (consumes remote runner minutes) | **Zero (local CPU cycles, < 200 ms for static audit)** |
| **Developer Flow** | Context-switching while waiting for remote checks | **Seamless (immediate feedback within active turn)** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Wiring CSpell, `check_links.py`, and `verify_traceability.py --audit` into `.githooks/pre-commit` shifts verification left into the agent's inner loop. If an agent introduces a broken link, typo, or malformed requirement tag, `git commit` fails immediately in the terminal output. The agent can self-correct in the same turn. The `.githooks/pre-push` hook executes `mise run pre:push` (`lint` + `test:affected`), guaranteeing that affected component and requirement tests pass before code leaves the workstation.

---

## 5. Target Architecture and Workflow

![RFL CI Execution Planner and Change Impact Analyser](images/ep5_ci_planner_c4.svg)

*Figure 5.1: C4 component architecture of the CI Execution Planner and Change Impact Analyser. Source script co-located at [`docs/eps/images/ep5_ci_planner_c4.d2`](images/ep5_ci_planner_c4.d2).*

### 5.1 Core Subsystem Architecture and Component Responsibilities
As illustrated in Figure 5.1, the traceability and execution planning pipeline divides responsibilities across five distinct architectural stages:

1. **Dependency Graph Builder:**
   Parses C++ `#include` statements across `src/` and `tests/` using an in-memory resolution cache to construct a bidirectional dependency graph without redundant filesystem access. Pairs `.cpp` implementation files with their corresponding `.hpp` headers. Traverses `src/python_bindings/bindings.cpp` to map exposed C++ types to `rfl.<Symbol>` exports. Parses Python test ASTs to map library imports and symbol calls back to underlying C++ components.

2. **Component and Test Resolver (Stage 1):**
   Ingests changed files from `git diff` against the merge base (or uncommitted working changes). Conducts a reverse breadth-first search (BFS) through the dependency graph to find all affected C++ compilation units and test executables. Traces Python test dependencies and flags build-wide modifications (such as `CMakeLists.txt` or `pyproject.toml`) that require full-suite execution.

3. **Requirement Impact Analyser (Stage 2):**
   Inspects the affected tests identified in Stage 1 to extract associated `[REQ-...]` tags. Scans `git diff` for modified requirement rows in `docs/eps/ep-*.md` specification tables. Parses `Implements: REQ-...` commit trailers to capture explicitly declared targets. Merges these inputs into a consolidated list of impacted requirements.

4. **Test Closure and Execution Planner (Stage 3):**
   Queries the traceability database to find all tests that verify the impacted requirements (requirement test closure). Merges directly affected tests with requirement-closure tests into a minimal, sufficient execution set. Emits structured execution plans (`plan.json`), regex filters for CTest (`ctest -R <filter>`), Python test target paths, and CI runner step outputs.

5. **Governance and Verification Gates:**
   Strictly separates static requirement mapping audits from dynamic execution evaluation. In static audit mode, the engine verifies that all claimed or impacted requirements exist in the specification and possess corresponding test cases in code. In dynamic verification mode, the engine inspects JUnit XML reports from CTest and Pytest to guarantee that tests passed without regression. Inverted control is strictly prohibited: the verification engine never executes test runners inline.

### 5.2 The Contributor Lifecycle

```
1. SPECIFY (Human Architect)
   - Define or refine requirements in docs/eps/: REQ-<SYS>-<NUM>
   - Specify physical invariant and numerical tolerance

2. IMPLEMENT AND TAG (Contributor / Agent)
   - Implement C++ / Python feature and unit tests (zero tags in src/)
   - Tag test cases with // [REQ-...] or @pytest.mark.req("REQ-...")
   - Local pre-commit hook verifies spelling, links, and traceability (< 1s)

3. VERIFY LOCALLY (Pre-Push Gate)
   - Local pre-push hook runs `mise run pre:push` (lint + affected component/requirement tests)

4. PLAN, QUALIFY AND MERGE (3-Tier CI and Human Review)
   - Tier 1 (`plan`): `mise run ci:plan` outputs `plan.json` (affected tests and impacted REQs)
   - Tier 2 (Execution): Runs `mise run ci:docs`, targeted CTest (`mise run ci:test:cpp`), and `mise run ci:verify:pr`
   - Tier 3 (`ci-gate`): `mise run ci:gate` verifies all `required_checks` in `plan.json` succeeded
```

### 5.3 CLI Interface (`mise` Contract Layer)

All local hooks and GitHub Actions workflows invoke `mise` targets rather than calling underlying scripts directly:

```bash
# Inspect affected components, test suites, and impacted requirements
mise run affected

# Generate deterministic CI execution plan as JSON
mise run ci:plan

# Run C++ and Python unit tests strictly on affected components and impacted requirements
mise run test:affected

# Execute C++ unit tests with optional regex filter and parallel execution
mise run test:cpp -- -R "^(DiracOperatorTests)\."

# Execute C++ test suite in CI runners with regex filter and JUnit XML output
mise run ci:test:cpp -- -C Release -R "^(CliffordTests)\." --junit build/junit.xml

# Audit repository for duplicate IDs, orphaned tags, and missing test coverage
mise run lint:traceability

# Print complete requirements-to-tests traceability matrix
mise run traceability:matrix

# Query detailed status, invariant, and mapped tests for a requirement
mise run traceability:query -- REQ-ARCH-001

# Validate pull request branch against origin/main (impacted REQs, test diffs, JUnit XML)
mise run ci:verify:pr
```

---

## 6. Verification Matrix and Quality Gates

| Verification Gate | Command / Target | Target Invariant |
| :--- | :--- | :--- |
| **Pre-Commit Hook** | `mise run pre:commit` | CSpell, markdown link checks, and traceability audit complete in ~1.5s. |
| **Pre-Push Hook** | `mise run pre:push` | Executes `lint` and `test:affected` before push. |
| **CI Execution Planner** | `mise run ci:plan` | Computes affected C++ and Python tests and impacted `REQ-` IDs in < 200 ms (`plan.json`). |
| **Static ID Audit** | `mise run lint:traceability` | All `REQ-` IDs are globally unique across EPs; zero orphaned test tags. |
| **Tool Unit Tests** | `mise run test:scripts` | Unit tests verify EP tables, test extraction, `pybind11` AST bridge, and `ci:gate`. |
| **CI Smart Gatekeeper** | `mise run ci:gate` | Enforces all `required_checks` from `plan.json` and renders impact tables in `$GITHUB_STEP_SUMMARY`. |

