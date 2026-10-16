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

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Standardise requirement namespacing and resolve ID collisions across EP-1 to EP-4 | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Complete |
| **Phase 2** | Implement zero-dependency traceability engine (`scripts/verify_traceability.py`) and `mise` tasks | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Complete |
| **Phase 3** | Integrate shift-left local Git hooks and CI Pull Request Step Summaries | `v0.3.0` | [#84](https://github.com/pauldruce/RFL/pull/84) | ✅ Complete |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Agentic Context
Software development in `RFL` pairs a human physicist-architect with autonomous AI coding agents.
Without clear structural guardrails, agentic workflows exhibit three recurring failure modes:

1. **Requirement Collision & Ambiguity:** Un-namespaced identifiers (such as `REQ-001` appearing in EP-1, EP-3, and EP-4) prevent deterministic lookup across specifications and test suites.
2. **Traceability Drift:** Scientific invariants defined in Enhancement Proposals (EPs) drift from the C++ and Python test suites because no automated tool links specifications to executable tests.
3. **Slow Feedback Loops:** Relying solely on remote CI pipelines adds a 2-to-5-minute delay per iteration. Agents finish their active turn before discovering broken links, spelling errors, or missing test tags.

### 2.2 Goals
* **Deterministic Verification:** Provide a zero-dependency Python tool (`scripts/verify_traceability.py`) that maps Markdown requirements (`REQ-<SYS>-<NUM>`) to GoogleTest, Catch2, and Pytest test cases.
* **Domain-Aware Auditing:** Distinguish between unit-testable code/physics invariants (`ARCH`, `PHY`, `TRC`) and CI/packaging workflow rules (`CI`, `PKG`) to eliminate false-positive audit warnings.
* **Shift-Left Agentic Feedback:** Execute fast static checks in `.githooks/pre-commit` (< 1 second) and full test suites in `.githooks/pre-push` so agents self-correct within their active turn.
* **Low-Friction PR Governance:** Publish a clear Markdown traceability table and modified-test diff alert to `$GITHUB_STEP_SUMMARY` without blocking internal refactoring PRs.
* **Compliance with EP-3:** Ensure all traceability checks execute in under 200 ms, preserving the sub-120-second CI budget (`REQ-CI-001`).

### 2.3 Non-Goals
* Do not adopt heavyweight external ALM frameworks (Sphinx-Needs, OpenFastTrace, ReqIF, or Doorstop).
* Do not enforce rigid multi-commit Git history choreography (such as mandatory two-commit Red-Green git time-travel) that breaks during rebasing or interactive C++ debugging.
* Do not require cryptographic SHA-256 manifest files or special GitHub override labels to edit or add tests.
* Do not track auto-generated traceability matrix files in Git where they cause merge conflicts across branches.

---

## 3. Requirements & Operational Invariants

| Requirement ID | Requirement Summary | Lifecycle Status | Operational & Architectural Invariant |
| :--- | :--- | :--- | :--- |
| **REQ-TRC-001** | **Global Requirement Namespacing** | Phase 1 (Implemented) | Requirement IDs must conform to regex `^REQ-[A-Z]{2,4}-\d{3}[a-z]?$`. Subsystem codes: `ARCH` (Architecture), `PKG` (Packaging), `CI` (CI/CD), `PHY` (Physics), `TRC` (Traceability). |
| **REQ-TRC-002** | **Lifecycle & Domain Awareness** | Phase 1 (Implemented) | Requirement tables in Section 3 of all EPs must include four columns: `Requirement ID`, `Requirement Summary`, `Lifecycle Status`, and `Operational & Physical Invariant`. Test coverage audits evaluate unit-testable requirements (`ARCH`, `PHY`, `TRC`) marked `Implemented`. |
| **REQ-TRC-003** | **Zero-Dependency Traceability Engine** | Phase 2 (Implemented) | `scripts/verify_traceability.py` must use standard Python $\ge 3.10$ libraries only. Exit codes: 0 (pass), 1 (validation error), 2 (CLI usage error). Runs in under 200 ms. |
| **REQ-TRC-004** | **Multi-Framework Tag Support** | Phase 2 (Implemented) | The engine must extract tags from GoogleTest comments (`// [REQ-...]`), Catch2 tags (`[REQ-...]`), and Pytest test functions (`def test_*` via `@pytest.mark.req`, docstrings, or comments), ignoring non-test helper methods. |
| **REQ-TRC-005** | **PR Step Summary & Test-Diff Visibility** | Phase 3 (Implemented) | PRs supplying `Implements: REQ-...` trailers are verified against static test mappings and CTest/Pytest JUnit XML outputs. Any modifications to existing test files are highlighted in `$GITHUB_STEP_SUMMARY`. |
| **REQ-TRC-006** | **Shift-Left Local Git Hooks** | Phase 3 (Implemented) | `.githooks/pre-commit` executes CSpell, link validation, and `verify_traceability.py --audit` on staged changes in under 1 second. `.githooks/pre-push` runs `mise run pre-push`. |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: Custom Python Standard Library Tooling vs Heavyweight ALM Frameworks

| Criteria | Option A: External Tools (Sphinx-Needs / Doorstop) | Option B: Custom Standard Library Python Script (Selected) |
| :--- | :--- | :--- |
| **Dependencies** | Heavy (Requires Sphinx build or external Python packages) | **Zero (Standard Python 3 only, native to RFL)** |
| **Execution Time** | 5–30 seconds in CI | **< 200 ms (instant local feedback)** |
| **Markdown Fit** | Poor (Requires YAML files or rST directives) | **Native (Parses GitHub Flavoured Markdown tables directly)** |
| **Maintenance** | Vulnerable to upstream deprecation | **Complete repository ownership alongside existing scripts** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* RFL follows lean development principles. Adding external documentation build chains violates `REQ-CI-001` and `REQ-CI-005`. A single-file Python script in `scripts/` mirrors existing repository utilities (`scripts/check_links.py`, `scripts/verify_changelog.py`).

---

### 4.2 ADR-2: JUnit Execution Proof & Diff Visibility vs Git-Topology & SHA-256 Ceremonies

| Criteria | Option A: Two-Commit Git Choreography & SHA-256 Manifests | Option B: JUnit XML Verification & PR Test-Diff Visibility (Selected) |
| :--- | :--- | :--- |
| **C++ Iteration Ergonomics** | ❌ Poor (Fights interactive debugging, rebasing, and fixup commits) | ✅ **Optimal (Standard atomic commits across headers, source, and tests)** |
| **Test Addition Friction** | ❌ High (Adding a test or fixing a comment breaks SHA-256 seals) | ✅ **Zero (New tests require no manifest regeneration or PR labels)** |
| **Tamper Visibility** | Cryptographic hash failure on file change | **Explicit alert in `$GITHUB_STEP_SUMMARY` listing modified test files** |
| **Execution Proof** | Fragile git checkout gymnastics in CI | **Direct verification of CTest and Pytest JUnit XML results** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Enforcing multi-commit Red-Green git topology or SHA-256 test manifests introduces high rigidity and encourages agents to write compliance theatre (such as dummy stubs or trivial arithmetic tests). In RFL, a human architect reviews every pull request. Verifying that claimed requirements (`Implements: REQ-...`) execute and pass in JUnit XML—while surfacing any modified existing test files directly in `$GITHUB_STEP_SUMMARY`—provides complete assurance with zero ceremonial overhead.

---

### 4.3 ADR-3: Informative PR Summaries vs Hard Blocking on Untagged PRs

| Criteria | Option A: Hard CI Blocker on Untagged PRs | Option B: Dynamic Informational Step Summary (Selected) |
| :--- | :--- | :--- |
| **Refactoring Friction** | High (Forces fake `REQ-` tags on routine cleanups) | **Zero (Untagged PRs pass as internal maintenance)** |
| **Visibility** | Obscure CI console log | **Clean Markdown table posted to `$GITHUB_STEP_SUMMARY`** |
| **Gating Boundary** | All PRs blocked unless tagged | **Hard gate only when `Implements: REQ-...` is claimed; informational otherwise** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Blocking commits that lack a formal requirement causes requirement inflation. If a PR supplies `Implements: REQ-...`, CI verifies that the requirement exists, maps to at least one test, and passes in JUnit XML. Routine refactoring PRs pass without friction.

---

### 4.4 ADR-4: Shift-Left Local Git Hooks vs Remote-Only CI Gating

| Criteria | Option A: Remote CI Only | Option B: Shift-Left Dual Loop (Local Hooks + Remote CI) (Selected) |
| :--- | :--- | :--- |
| **Agent Feedback Latency** | 2–5 minutes per push | **< 1 second (instant stderr self-correction in local turn)** |
| **PR Commit Hygiene** | Poor (noisy retry commits) | **Optimal (clean commits verified before push)** |
| **Compute Consumption** | High (consumes remote runner minutes) | **Zero (local CPU cycles, < 200 ms for static audit)** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Wiring CSpell, `check_links.py`, and `verify_traceability.py --audit` into `.githooks/pre-commit` shifts verification left into the agent's inner loop. If an agent introduces a broken link, typo, or malformed requirement tag, `git commit` fails immediately in the terminal output, allowing the agent to self-correct in the same turn. `.githooks/pre-push` executes `mise run pre-push` to guarantee all unit tests pass before code leaves the workstation.

---

## 5. Target Architecture & Workflow

### 5.1 The Contributor Lifecycle

```
1. SPECIFY (Human Architect)
   - Define or refine requirement in docs/eps/: REQ-<SYS>-<NUM>
   - Specify physical invariant and numerical tolerance

2. IMPLEMENT & TAG (Contributor / Agent)
   - Implement C++ / Python feature and unit tests
   - Tag test cases with // [REQ-<SYS>-<NUM>] or @pytest.mark.req("REQ-...")
   - Local pre-commit hook verifies spelling, links, and traceability (< 1s)

3. VERIFY & PUSH (Local Pre-Push Gate)
   - Include "Implements: REQ-<SYS>-<NUM>" in commit message or PR body
   - Local pre-push hook runs `mise run pre-push` (linting + CTest + Pytest)

4. QUALIFY & MERGE (CI & Human Review)
   - CI validates requirement mappings and JUnit XML execution
   - CI posts traceability table and modified-test alerts to $GITHUB_STEP_SUMMARY
   - Human architect reviews diff and merges via Squash-and-Merge
```

### 5.2 CLI Interface (`scripts/verify_traceability.py`)

```bash
# Audit repository for duplicate IDs, orphaned tags, and testable coverage
python3 scripts/verify_traceability.py --audit

# Query requirement details, invariant, and mapped test cases
python3 scripts/verify_traceability.py --query REQ-ARCH-001

# Validate current PR branch against origin/main (trailers, test diffs, JUnit XML)
python3 scripts/verify_traceability.py --check-pr --base origin/main [--verify-xml build/junit.xml]

# Generate Markdown traceability matrix on demand (stdout or file)
python3 scripts/verify_traceability.py --generate-matrix [docs/Traceability_Matrix.md]

# Dump machine-readable JSON representation for AI agents
python3 scripts/verify_traceability.py --dump-json
```

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Target | Target Invariant |
| :--- | :--- | :--- |
| **Pre-Commit Hook** | `.githooks/pre-commit` | CSpell, markdown link checks, and `verify_traceability.py --audit` complete in < 1s. |
| **Pre-Push Hook** | `.githooks/pre-push` | Runs `mise run pre-push` (all lint checks, CTest, and Pytest) before pushing. |
| **Static ID Audit** | `mise run lint-traceability` | All `REQ-` IDs are globally unique across EPs; zero orphaned test tags. |
| **Tool Unit Tests** | `pytest scripts/tests/test_verify_traceability.py` | Unit test coverage over Markdown table, GTest, Catch2, Pytest AST, and JUnit XML parsing. |
| **CI PR Summary** | `.github/workflows/ci.yml` | Renders claimed requirement status and modified test file alerts in `$GITHUB_STEP_SUMMARY`. |
