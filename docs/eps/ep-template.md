# EP-<N>: [Feature / Enhancement Title]

* **Title:** [Feature / Enhancement Title]
* **Author:** [Author Name]
* **Status:** [Draft / In Discussion / Accepted (In Progress) / Completed]
* **Target Versions:** [e.g. RFL v0.1.0 (Phase 1), v0.3.0 (Phase 2)]
* **Date:** YYYY-MM-DD

---

## 1. Multi-Phase Implementation Tracker

This proposal spans multiple releases.
The table below tracks the status of each implementation phase:

| Phase | Scope & Deliverables | Target Version | PR / Issue | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Phase 1** | [Deliverables summary] | `v0.1.0` | [PR #X](...) (Closes [#Y](...)) | ⏳ Scheduled |
| **Phase 2** | [Deliverables summary] | `v0.3.0` | [#Z](...) | 💡 Planned |

---

## 2. Motivation, Goals & Non-Goals

### 2.1 Problem Statement & Research Context
[Describe the physics, mathematics, or computational bottleneck being addressed.]

### 2.2 Goals
* **Goal 1:** [Specific research or technical capability.]
* **Goal 2:** [Specific scientific invariant or optimisation.]

### 2.3 Non-Goals
* [Explicitly declare what is out of scope for this proposal.]

---

## 3. Research Workflows & Scientific Requirements

### 3.1 Core Research Scenarios
1. **Scenario 1 (Interactive Exploration in Python):** [Describe researcher interactive workflow.]
2. **Scenario 2 (Automated Batch Sampling):** [Describe batch simulation and statistics workflow.]

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Physical & Mathematical Invariant |
| :--- | :--- | :--- |
| **REQ-001** | **[Requirement Name]** | [Mathematical or physical invariant to preserve.] |
| **REQ-002** | **[Requirement Name]** | [Computational precision or tolerance limit.] |

---

## 4. Architecture Decision Records (ADRs) & Trade-offs

### 4.1 ADR-1: [Subsystem / Architectural Decision]

| Criteria | Option A: [Option Name] (Selected) | Option B: [Option Name] |
| :--- | :--- | :--- |
| **Performance** | High | Low |
| **API Ergonomics** | Clean value type | Complex heap wrapper |
| **Decision** | **Selected (Option A)** | Rejected |

*Rationale:* [Clear justification for the selected approach using ASD-STE100 controlled English.]

---

## 5. Target Architecture & Component Contracts

### 5.1 Subsystem Architecture
[Define C++ class declarations, signatures, value types, and file organisation.]

### 5.2 Python Interoperability
[Define pybind11 bindings, zero-copy buffer views, and NumPy integration.]

---

## 6. Verification Matrix & Quality Gates

| Verification Gate | Command / Test Description | Target Invariant |
| :--- | :--- | :--- |
| **Unit Tests** | `ctest --test-dir build` | Verifies individual component contracts. |
| **Integration Tests** | `pytest src/RFL/python_bindings/tests` | Verifies end-to-end Python bindings and buffer safety. |

---

## 7. Phased Delivery Plan

### Phase 1: [Phase 1 Title]
* **Target Version:** `vX.Y.Z`
* **GitHub Issue:** [Issue #A](...)
* **Tasks:**
  1. [Task 1]
  2. [Task 2]

### Phase 2: [Phase 2 Title]
* **Target Version:** `vX.Y.Z`
* **GitHub Issue:** [Issue #B](...)
* **Tasks:**
  1. [Task 1]
  2. [Task 2]
