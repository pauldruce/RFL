# EP-5: Specification-Driven Agentic Workflow and Traceability Infrastructure

* **Title:** Specification-Driven Agentic Workflow and Traceability Infrastructure
* **Author:** Paul Druce
* **Status:** Accepted
* **Target Versions:** RFL v0.3.0 (Phase 1 and Phase 2), v0.4.0 (Phase 3 and Phase 4)
* **Date:** 2026-10-14

---

## 1. Multi-Phase Implementation Tracker

This proposal establishes the workflow and tooling infrastructure for autonomous AI agent development and human qualification.
The table below tracks the status of each implementation phase:

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | Standardise requirement namespacing and resolve ID collisions across EP-1 to EP-4 | `v0.3.0` | — | ⏳ Scheduled |
| **Phase 2** | Implement zero-dependency traceability engine (`scripts/verify_traceability.py`) and `mise` tasks | `v0.3.0` | — | 💡 Planned |
| **Phase 3** | Implement physical invariant test sealing (`tests/qualification/.manifest` SHA-256 integrity) | `v0.4.0` | — | 💡 Planned |
| **Phase 4** | Integrate GitHub Actions CI PR Step Summaries and two-commit Red-Green verification | `v0.4.0` | — | 💡 Planned |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Agentic Context
Software development in `RFL` increasingly uses autonomous AI coding agents as primary implementers, with humans acting as architects, specifiers, and qualification reviewers. 

This paradigm introduces specific operational risks:
1. **Architectural Erosion & Scope Creep:** Agents take the path of least resistance. They skip layer boundaries, introduce circular dependencies, or modify unrelated files across the repository.
2. **Tautological Tests & False Confidence:** When instructed to follow Test-Driven Development (TDD), agents often implement code first and write trivial tests that pass immediately without verifying the underlying physics.
3. **The "Passing Fix" & Compilation Deadlocks:** Rigid multi-PR path gates forbid implementation PRs from modifying test files. This causes fatal C++ compilation errors during API refactoring and blocks bug fixes when a corrected test passes against `main`.
4. **Traceability Vacuum:** Scientific requirements in Enhancement Proposals (EPs) drift from code. There is no automated mechanism to prove that an implemented feature satisfies its mathematical specification.

### 2.2 Goals
* **Deterministic Verification:** Provide a zero-dependency Python tool (`scripts/verify_traceability.py`) that maps Markdown requirements (`REQ-<SYS>-<NUM>`) to GoogleTest, Catch2, and Pytest suites.
* **Atomic Two-Commit Verification:** Verify TDD compliance (Red on stubs, Green on implementation) within a single atomic PR without compilation deadlocks.
* **Qualification Invariant Sealing:** Prevent silent degradation of numerical tolerances by sealing acceptance test files with cryptographic checksums.
* **Informative CI Governance:** Generate automated Markdown compliance matrices in GitHub Pull Request Step Summaries without blocking internal refactoring tasks.
* **Compliance with EP-3:** Ensure all traceability and verification checks execute in under 2 seconds, preserving the sub-120-second CI budget (`REQ-CI-001`).

### 2.3 Non-Goals
* Do not adopt heavyweight external ALM or traceability frameworks (Sphinx-Needs, OpenFastTrace, ReqIF, or Doorstop).
* Do not forbid atomic C++ commits across `include/`, `src/`, and `tests/`.
* Do not run exhaustive mutation testing during pull request gates.

---

## 3. Requirements & Operational Invariants

| Requirement ID | Requirement Summary | Lifecycle Status | Operational & Architectural Invariant |
| :--- | :--- | :--- | :--- |
| **REQ-TRC-001** | **Global Requirement Namespacing** | Phase 1 (Planned) | Requirement IDs must conform to regex `^REQ-[A-Z]{2,4}-\d{3}$`. Subsystem codes: `ARCH` (Architecture), `PKG` (Packaging), `CI` (CI/CD), `PHY` (Physics), `TRC` (Traceability). |
| **REQ-TRC-002** | **Lifecycle Phase Awareness** | Phase 1 (Planned) | Requirement tables in Section 3 of all EPs must include four columns: `Requirement ID`, `Requirement Summary`, `Lifecycle Status`, and `Operational / Physical Invariant`. The engine must evaluate only requirements marked `Phase <N> (Implemented)`. |
| **REQ-TRC-003** | **Zero-Dependency Traceability Engine** | Phase 2 (Planned) | `scripts/verify_traceability.py` must use standard Python $\ge 3.10$ libraries only. Exit codes: 0 (pass), 1 (traceability error), 2 (CLI usage/parser error). Runs using `mise` in under 200 ms. |
| **REQ-TRC-004** | **Multi-Framework Tag Support** | Phase 2 (Planned) | The tool must extract tags from GoogleTest comments (`// [REQ-...]`) preceding test macros, Catch2 bracketed tags (`[REQ-...]`), and Pytest AST decorators (`@pytest.mark.req`). |
| **REQ-TRC-005** | **Conventional PR Metadata & Dynamic Verification** | Phase 2 (Planned) | Feature PRs claiming requirement delivery must supply Git trailers (`Implements: REQ-...`). The tool maps static source tags to test symbols, and confirms in JUnit XML that all target tests executed with zero failures. Catch2 v3 and Pytest omit `status="run"`; verification asserts the test element is present and contains no `<skipped>`, `<error>`, or `<failure>` tags. |
| **REQ-TRC-006** | **Qualification Test Sealing** | Phase 3 (Planned) | Core scientific qualification tests in `tests/qualification/` must match SHA-256 hashes in `.manifest` using normalised LF line endings. Tolerance changes require `CODEOWNERS` approval and the `qualification-override` PR label. |
| **REQ-TRC-007** | **Atomic Red-Green Protocol** | Phase 4 (Planned) | Feature PRs introducing new components must provide a Red Anchor and Green Anchor (`HEAD`). The Red Anchor is resolved as the parent of the first `feat:` commit (`git rev-list --reverse --parents -E --grep='^feat(\([^)]+\))?:' origin/main..HEAD`). Stubs return degenerate values (`NaN` or zeros); at the Red Anchor, all tests mapped to the PR's `Implements: REQ-...` trailer(s) must produce `<failure>`, all pre-existing unrelated tests must pass, and the suite must produce zero `<error>` tags. CI runs an in-place checkout in the same build directory on Ubuntu. Repositories enforce Squash-and-Merge to protect `git bisect` on `main`. |
| **REQ-TRC-008** | **Fortran Alignment & Zero-Copy Invariant** | Phase 2 (Planned) | Python buffer bindings over Armadillo matrices must verify column-major / Fortran alignment (`.flags['F_CONTIGUOUS']`) and pointer equality. This prevents silent matrix transposition and complex conjugation ($H \to H^*$) of Hermitian operators. |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: Custom Python Standard Library Tooling vs Heavyweight ALM Frameworks

| Criteria | Option A: External Tools (Sphinx-Needs / OpenFastTrace) | Option B: Custom Standard Library Python Script (Selected) |
| :--- | :--- | :--- |
| **Dependencies** | Heavy (Requires Sphinx documentation build or Java Runtime Environment) | **Zero (Standard Python 3 only, native to RFL)** |
| **Execution Time** | 10–30 seconds in CI | **< 200 ms (instant local feedback)** |
| **Markdown Fit** | Poor (Requires external specification files, YAML fragmentation, or rST directives) | **Native (Parses GitHub Flavoured Markdown tables directly)** |
| **Maintenance** | Vulnerable to upstream tool deprecation and version mismatches | **Complete repository ownership alongside existing verification tools** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* RFL adheres strictly to lean development principles. Adding Java or Sphinx runtimes violates `REQ-CI-001` and `REQ-CI-005`. A single-file Python script in `scripts/` mirrors existing successful verification scripts (`scripts/check_d2_images.py`, `scripts/verify_changelog.py`).

---

### 4.2 ADR-2: Atomic Red-Green Protocol vs Rigid Multi-PR Path Gating

| Criteria | Option A: Rigid Split PRs (`spec:` -> `test:` -> `feat:`) | Option B: Atomic PR with In-Place Red-Green Protocol (Selected) |
| :--- | :--- | :--- |
| **C++ Refactoring** | ❌ Fails (Breaks atomic AST compilation during signature changes) | ✅ **Optimal (Interface, tests, and call sites update atomically)** |
| **Bug Fix Tolerance** | ❌ Fails (Passing fixes cannot merge if tests must fail against `main`) | ✅ **Optimal (Tagged as `fix:` or tolerance recalibration)** |
| **TDD Verification** | Requires multi-PR reviews and manual dependency coordination | ✅ **Automated (Red Anchor checked for assertion `<failure>`; Green Anchor checked for pass)** |
| **CI Build Time** | 6–10 minutes (Double checkout and recompilation across 3 OSs) | ✅ **< 10 seconds on Ubuntu runner (In-place checkout with 100% ccache hit and incremental Ninja)** |
| **Template Support** | ❌ Fails (Inline expression templates cannot be modified in Commit 2) | ✅ **Supported (Semantic AST check: inline bodies update if public signatures match)** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Red-team analysis proved that rigid multi-PR path gates introduce fatal deadlocks in statically compiled C++. Enforcing an atomic PR with a Red Anchor and Green Anchor (`HEAD`) guarantees TDD verification without breaking atomic compilation. 

The Red Anchor is resolved deterministically as the **immediate parent of the first `feat:` commit** (`git rev-list --reverse --parents -E --grep='^feat(\([^)]+\))?:' origin/main..HEAD | head -n 1 | awk '{print $2}'`). This represents the exact state of the repository before implementation began. It permits contributors to push multiple iterative `test:` commits (e.g. fixing syntax errors or missing headers) without failing CI. Alternate strategies were evaluated and rejected:
* *Resolution 1 (Last `test:` commit):* Rejected due to **chronological inversion**. If an agent pushes an edge-case test fix after `feat:`, the last `test:` commit already contains the implementation, causing the Red check to pass falsely and fail CI.
* *Resolution 3 (AST body-stripping via tree-sitter):* Rejected because tree-sitter lacks semantic C++ type resolution. Replacing Armadillo function bodies with `return {};` generates $0 \times 0$ matrices, triggering out-of-bounds indexing and dimension mismatch exceptions (`<error>`) instead of assertion failures (`<failure>`). It also breaks on reference return types (`const arma::cx_mat&`), `auto` return deduction, and constructor initialiser lists.

Stubs must return **degenerate values** (such as `arma::zeros(...)` or `NaN`), **never throwing uncaught exceptions**. In JUnit XML, an uncaught exception produces `<error>`, which aborts test execution before assertions run; degenerate values force the test to execute and fail a genuine mathematical assertion (`<failure>`). CI runs an **in-place checkout** in the single build directory on Ubuntu: `git checkout <red_sha> && cmake --build build && ctest`, followed by `git checkout <green_sha> && cmake --build build && ctest`. This preserves identical working directory hashes for `ccache` and Ninja logs, recompiling only modified translation units in under 2 seconds. Repositories enforce Squash-and-Merge so `main` never receives a failing commit.

---

### 4.3 ADR-3: Informative PR Summaries vs Blocking Path Gates

| Criteria | Option A: Hard CI Blocker on Untagged PRs | Option B: Dynamic Informational Step Summary (Selected) |
| :--- | :--- | :--- |
| **Refactoring Friction** | High (Forces developers to attach fake `REQ-` tags to routine cleanups) | **Zero (Untagged PRs marked as internal refactoring)** |
| **Visibility** | Obscure CI console log | **Clean Markdown table posted to** `$GITHUB_STEP_SUMMARY` |
| **Execution Proof** | Static string matching only (Vulnerable to disabled/skipped tests) | **Dynamic parsing of CTest and Pytest JUnit XML test results** |
| **Gating Boundary** | All PRs blocked unless tagged | **Hard gate on claimed requirements and sealed test modifications; non-blocking on refactoring** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Blocking every commit that lacks a formal requirement causes requirement inflation and developer circumvention. Dynamic verification of declared requirements using JUnit XML gives complete transparency. If a PR supplies `Implements: REQ-...`, CI blocks if the test is missing or failing; untagged refactoring PRs remain informational.

---

### 4.4 ADR-4: Shift-Left Local Git Hooks vs Remote-Only CI Gating

| Criteria | Option A: Remote CI Only | Option B: Shift-Left Dual Loop (Local Hooks + Remote CI) (Selected) |
| :--- | :--- | :--- |
| **Agent Feedback Latency** | 2–5 minutes per push | **< 1 second (instant stderr self-correction in local turn)** |
| **PR Commit Hygiene** | Poor (multiple noisy retry commits) | **Optimal (clean, atomic commits on branch)** |
| **Context Retention** | Low (agent ends turn prematurely) | **High (agent self-corrects in the active execution turn)** |
| **Compute Consumption** | High (consumes remote runner minutes) | **Zero (local CPU cycles, < 200 ms)** |
| **Integrity Assurance** | Authoritative (clean-room multi-OS) | **Authoritative (CI mirrors identical checks on clean virtual machines)** |
| **Decision** | Rejected | **Selected (Option B)** |

*Rationale:* Autonomous AI coding agents optimise for rapid task completion. When verification happens solely in remote CI pipelines, feedback latency averages 2 to 5 minutes. The agent concludes its execution turn falsely, leaving broken links or failing tests for the user to report.

Wiring `verify_traceability.py --audit` and CSpell into `.githooks/pre-commit` shifts verification left into the agent's inner loop. The `git commit` command fails immediately in the terminal tool output. The agent observes the error in stderr and corrects the issue during the active turn without human intervention. The `.githooks/pre-push` hook executes `mise run pre-push`, preventing broken builds from leaving the workstation. Remote CI remains essential as the authoritative clean-room gatekeeper across Ubuntu, macOS, and Windows.

---

## 5. Target Architecture & The Specification-Driven Development Workflow

### 5.1 The 5-Stage Contributor Lifecycle

The workflow operates as an identity-agnostic pipeline. It applies equally to human developers and autonomous AI agents:

```
Stage 1: SPECIFY (Human Architect)
  │  - Author requirement in EP: REQ-<SYS>-<NUM>
  │  - Define invariants and numerical tolerances
  ▼
Stage 2: SCAFFOLD (Contributor / Agentic Loop)
  │  - Declare public headers in include/rfl/ with degenerate stubs (zeros/NaN)
  │  - Register tagged tests (// [REQ-...] or @pytest.mark.req)
  │  - Local pre-commit hook verifies requirement tags and spelling (< 1s)
  ▼
Stage 3: IMPLEMENT VIA RED-GREEN ANCHORS (Contributor / Agentic Loop)
  │  - Red Anchor (test: commit): Compiles, links, and produces assertion <failure> (Red)
  │  - Green Anchor (HEAD): Replaces stubs with math kernel (Green)
  │  - Append Git trailer: "Implements: REQ-<SYS>-<NUM>"
  │  - Local pre-push hook executes mise run pre-push (lint + unit tests)
  ▼
Stage 4: QUALIFY (Automated CI Pipeline)
  │  - Verify Git trailers against docs/eps/
  │  - In-place checkout on Ubuntu: verify Red <failure> on Red Anchor, Green on HEAD
  │  - Parse JUnit XML to confirm test execution (reject skipped tests)
  │  - Post Markdown compliance matrix to $GITHUB_STEP_SUMMARY
  ▼
Stage 5: SIGN-OFF & MERGE (Human Architect)
     - Audit Step Summary compliance matrix
     - Merge pull request using Squash-and-Merge
```

#### Detailed Stage Rules:

1. **Stage 1 (Specify):** The human architect authors or refines the Enhancement Proposal in `docs/eps/`. The EP defines the mathematical requirements, tolerances, and ADRs.
2. **Stage 2 (Scaffold):** The contributor declares public headers in `include/rfl/` using regular value types. The contributor writes test cases tagged with requirement identifiers (`// [REQ-...]` or `@pytest.mark.req`). Method stubs return degenerate values (`arma::zeros(...)` or `NaN`), never throwing uncaught exceptions. The `.githooks/pre-commit` hook audits requirement ID syntax, CSpell dictionaries, and markdown links in under one second.
3. **Stage 3 (Implement via Red-Green Anchors):**
   * The contributor creates an atomic pull request containing at least one commit matching `^test:` (the Red Anchor) and finishing at `HEAD` (the Green Anchor). The Red Anchor is deterministically resolved as the immediate parent of the first `feat:` commit.
   * **Red Anchor (`test:`):** Adds interface stubs and the new test assertions. It must compile and link with zero errors. All tests mapped to the claimed requirements must fail on mathematical assertions (`<failure>`), all pre-existing tests must pass, and the suite must produce zero uncaught exceptions (`<error>`).
   * **Green Anchor (`HEAD`):** Replaces stubs with the implementation in `src/`. All tests turn green. Feature commits may modify `include/` for templates and inline helpers provided public AST signatures remain identical.
   * The commit message or PR description must include the Git trailer: `Implements: REQ-...`.
   * The `.githooks/pre-push` hook executes `mise run pre-push` (incremental build, lint checks, and unit tests), preventing broken code from leaving the workstation.
4. **Stage 4 (Qualify):**
   * GitHub Actions runs `scripts/verify_traceability.py --check-pr`.
   * On the Ubuntu runner, CI executes an **in-place checkout**: it checks out the Red Anchor, compiles, and verifies that claimed tests produce `<failure>`, pre-existing tests pass, and zero `<error>` tags exist. It then checks out `HEAD` and recompiles incrementally using Ninja and `ccache` in under 2 seconds, verifying all tests pass.
   * CI checks JUnit XML reports to ensure claimed tests ran and passed (rejecting `<skipped>` tests).
   * The job posts the verification table to `$GITHUB_STEP_SUMMARY`.
5. **Stage 5 (Sign-Off & Merge):** The human reviewer audits the compliance matrix and architectural invariants. The maintainer merges the PR using Squash-and-Merge to maintain clean git history.

### 5.2 Component Structure

```
RFL Traceability & Governance Architecture
├── .githooks/
│   ├── pre-commit                                   (Shift-left hygiene: CSpell + verify_traceability)
│   └── pre-push                                     (Local smoke tests: mise pre-push)
├── docs/eps/
│   ├── ep-1-core-architecture-modernisation.md      (REQ-ARCH-001..007)
│   ├── ep-2-package-and-binary-distribution.md      (REQ-PKG-001..009)
│   ├── ep-3-resource-efficient-ci-pipeline.md       (REQ-CI-001..010)
│   ├── ep-4-computational-physics-verification-suite.md (REQ-PHY-001..011)
│   └── ep-5-specification-driven-agentic-workflow.md (REQ-TRC-001..008)
├── tests/
│   └── qualification/                               (Sealed acceptance tests)
│       └── .manifest                                (SHA-256 integrity hash)
├── src/
│   └── core/tests/                                  (Collocated unit tests)
├── scripts/
│   ├── verify_traceability.py                       (Core extraction & audit engine)
│   └── verify_tdd_commits.py                        (Incremental Red-Green runner)
└── .github/workflows/
    └── ci.yml                                       (Injects Step Summary)
```

### 5.3 Tooling CLI Interface (`scripts/verify_traceability.py`)

```bash
# Audit repository for duplicate IDs, orphaned tests, and coverage
python3 scripts/verify_traceability.py --audit

# Validate current PR branch against origin/main (Git trailers and test diff)
python3 scripts/verify_traceability.py --check-pr [--base origin/main]

# Query requirement summary and mapped test cases
python3 scripts/verify_traceability.py --query REQ-PHY-001

# Verify dynamic execution against multiple CTest and Pytest JUnit XML outputs
python3 scripts/verify_traceability.py --verify-xml build/junit-ctest.xml build/junit-pytest.xml

# Generate Markdown traceability matrix for documentation and Obsidian vault
python3 scripts/verify_traceability.py --generate-matrix docs/Traceability_Matrix.md

# Dump machine-readable JSON representation for AI agents
python3 scripts/verify_traceability.py --dump-json

# Verify integrity of sealed qualification tests using normalised line endings
python3 scripts/verify_traceability.py --verify-manifest tests/qualification/.manifest
```

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Target | Target Invariant |
| :--- | :--- | :--- |
| **Pre-Commit Shift-Left** | `.githooks/pre-commit` | CSpell spelling check and `verify_traceability.py --audit` execute on staged files in < 1 second. |
| **Pre-Push Smoke Gate** | `.githooks/pre-push` | Executes `mise run pre-push` (lint checks and fast unit tests) before remote push. |
| **Static ID Audit** | `mise run lint-traceability` | All `REQ-` IDs are globally unique; zero collisions across EPs. |
| **Unit Test Coverage** | `pytest scripts/tests/test_verify_traceability.py` | 100% test coverage over Markdown table, GTest, Catch2, and Pytest parsing. |
| **Execution Performance** | `time mise run lint-traceability` | Wall-clock execution time < 200 ms. |
| **CI Summary Generation** | GitHub Actions Workflow Step | Renders full traceability status in `$GITHUB_STEP_SUMMARY`. |
| **Integrity Check** | `python3 scripts/verify_traceability.py --verify-manifest tests/qualification/.manifest` | Qualification tests match sealed reference hashes using normalised LF line endings across OSs. |

---

## 7. Phased Delivery Plan

### Phase 1: Zero-Conflict Requirement Namespacing (Expand Phase)
* **Scope:** Resolve `REQ-001` collisions between EP-1 and EP-4. Establish standard prefixes:
  * `REQ-ARCH-xxx`: Core software architecture and value semantics (EP-1).
  * `REQ-PKG-xxx`: Packaging, ABI, and binary distribution (EP-2).
  * `REQ-CI-xxx`: Resource-efficient CI/CD infrastructure (EP-3, migrated from `REQ-ARC-xxx`).
  * `REQ-PHY-xxx`: Computational physics and numerical verification (EP-4).
  * `REQ-TRC-xxx`: Traceability and agentic workflow governance (EP-5).
* **Deliverables:**
  1. Update Markdown tables in `docs/eps/` with four-column schema and Phase status.
  2. **Zero-Conflict Aliasing:** Introduce in-memory alias mapping (`LEGACY_REQ_ALIASES: dict[tuple[str | None, str], str]`) inside `scripts/verify_traceability.py`. Colliding legacy IDs (`REQ-001` through `REQ-007`) are disambiguated by source file path prefix, mapping `(prefix, legacy_id)` tuple keys (such as `("src/core/tests/test_architecture", "REQ-001") -> "REQ-ARCH-001"` and `("src/core/tests/test_physics", "REQ-001") -> "REQ-PHY-001"`), with fallback to `(None, legacy_id)` for unambiguous legacy IDs. **Zero C++ or Python test files are modified on `main` during Phase 1**, guaranteeing zero merge conflicts for active in-flight feature branches (such as `feat/issue-14-barrett-glaser-action`).
  3. **Contract Phase:** Execute an automated refactoring script across `src/core/tests/` to update legacy test tags only after active feature branches have landed.

### Phase 2: Standalone Traceability Tooling
* **Scope:** Implement `scripts/verify_traceability.py` using Python standard library.
* **Deliverables:**
  1. Parser for EP Markdown requirement tables with LaTeX math pipe shielding.
  2. Extractor for GoogleTest, Catch2, and Pytest `ast.parse` annotations.
  3. Register `req(id)` custom test marker in `pyproject.toml` under `[tool.pytest.ini_options]`.
  4. Git trailer inspector (`Implements: REQ-...`).
  5. Unit test suite in `scripts/tests/test_verify_traceability.py`.
  6. Task registration in `mise.toml` (`lint-traceability`, `matrix-traceability`).

### Phase 3: Invariant Test Sealing Mechanism
* **Scope:** Establish `tests/qualification/` directory with `.manifest` SHA-256 seal.
* **Deliverables:**
  1. Enforce `eol=lf` for `tests/qualification/**` in `.gitattributes`.
  2. Add `CODEOWNERS` protection for `tests/qualification/`.
  3. Pre-commit and CI verification script ensuring scientific tolerances cannot be modified without maintainer cryptographic override.

### Phase 4: CI Step Summary & Two-Commit Protocol Integration
* **Scope:** Integrate dynamic JUnit XML inspection and TDD verification into `.github/workflows/ci.yml`.
* **Deliverables:**
  1. Automated Markdown table generation attached to GitHub PR summaries.
  2. Implementation of `scripts/verify_tdd_commits.py` using in-place checkout exclusively on the Ubuntu runner.
  3. Mandate Squash-and-Merge in repository branch settings.
