---
name: technical-writing
description: Use whenever writing, editing, or reviewing documentation, C++ docstrings, code comments, Enhancement Proposals (EPs), or research notes.
---

# Technical Writing & Controlled Language Guide

This guide defines the writing rules for RFL documentation, code comments, docstrings, and research notes. It adapts the international **ASD-STE100** specification with **British English spelling**.

> [!IMPORTANT]
> The single source of truth for approved terms and their definitions is **[docs/Glossary.md](docs/Glossary.md)**. Always check the glossary to maintain the "one word for one concept" rule.

---

## 1. Core Principles

Controlled writing ensures clarity, reduces ambiguity, and makes technical documentation easy to read for humans and AI agents.

| Pillar | Standard |
| :--- | :--- |
| **1. Short Sentences** | Maximum 20 words for instructions; maximum 25 words for descriptions |
| **2. Active Voice** | Direct verbs and imperatives |
| **3. Controlled Vocabulary** | One meaning per word |
| **4. British Spelling** | Standardise on British English (`-ise`, `-our`, `-re`, double `l`) |

---

## 2. The 8 Core Writing Rules

### Rule 1: Sentence Length Limits
* **Procedural steps & instructions:** Maximum **20 words** per sentence.
* **Descriptive text & explanations:** Maximum **25 words** per sentence.
* Keep one main idea per sentence.

### Rule 2: Use the Active Voice
* **Do not write (Passive):** *"The Dirac operator matrices are updated by the Metropolis algorithm."*
* **Write (Active):** *"The Metropolis algorithm updates the Dirac operator matrices."*

### Rule 3: Use Direct Imperative Verbs for Instructions
* **Do not write:** *"You should run the tests using ctest."*
* **Write:** *"Run the tests with `ctest`."*

### Rule 4: Avoid Complex Noun Clusters
* Limit noun sequences to a maximum of **3 nouns**.
* **Do not write:** *"Finite noncommutative geometry Dirac operator matrix element variation."* (7 nouns)
* **Write:** *"Variation of a matrix element in the finite NCG Dirac operator."*

### Rule 5: One Meaning per Word (No Ambiguous Conjunctions)
* Use **`because`** (not *as* or *since*) when giving a reason.
* Use **`after`** / **`when`** (not *as* or *since*) for time relationships.
* Use **`to`** (not *in order to*).
* Use **`must`** (not *shall*, *should*, or *ought to*) for mandatory requirements.

### Rule 6: Use Approved Verbs over Vague Words
| Vague / Unapproved | Approved Replacement | Example |
| :--- | :--- | :--- |
| *carry out / perform / conduct* | **do / run / execute / calculate** | *"Run the simulation."* (not *"Carry out the simulation."*) |
| *utilize / leverage* | **use** | *"Use Armadillo for matrix math."* |
| *terminate* | **stop / end / cancel** | *"Stop the iteration."* |
| *facilitate* | **help / enable / provide** | *"This method enables fast lookup."* |
| *via* | **with / through / using** | *"Sample using Metropolis."* |

### Rule 7: Avoid Continuous (-ing) Tenses in Procedures
* Prefer simple present or imperative over present continuous.
* **Do not write:** *"When running the sampler, it is recording eigenvalues."*
* **Write:** *"When the sampler runs, it records eigenvalues."*

---

## 3. British English Spelling Conventions

RFL strictly standardises on **British English**:

<!-- cspell:disable -->
| Feature | British Standard (Approved) | US Form (Avoid) |
| :--- | :--- | :--- |
| **-ise / -isation** | *initialise, randomise, optimise, diagonalisation, categorise* | *initialize, randomize, optimize, diagonalization* |
| **-our** | *behaviour, colour, neighbour* | *behavior, color, neighbor* |
| **-re** | *centre, metre, fibre* | *center, meter, fiber* |
| **Double 'l'** | *modelling, initialised, cancelled, travelling* | *modeling, initialized, canceled, traveling* |
| **-programme** | *programme* (for scientific initiatives), *program* (for computer code) | *program* |
| **-ence / -ense** | *licence* (noun), *license* (verb); *defence* | *license* (both); *defense* |
| **Specialised terms** | *gauge, analogue, catalogue* | *gage, analog, catalog* |
<!-- cspell:enable -->

---

## 4. Technical Writing Examples (Before vs. After)

### Example 1: Code Docstrings
* ❌ **Before:**
  ```cpp
  // This method is utilized in order to carry out the calculation of the trace of the fourth power
  // of the Dirac operator which is needed since we want to evaluate the Barrett-Glaser action.
  ```
* ✅ **After (ASD-STE100 + British):**
  ```cpp
  /**
   * Calculates the trace of the fourth power of the Dirac operator, Tr(D^4).
   * The Barrett-Glaser action uses this trace to evaluate energy.
   *
   * @return The trace value as a double.
   */
  ```

### Example 2: Architecture Documentation
* ❌ **Before:**
  ```markdown
  Since the previous implementation was utilizing dynamic polymorphism with unique_ptr wrappers,
  in order to achieve optimal computational performance we are refactoring the state into regular value types.
  ```
* ✅ **After (ASD-STE100 + British):**
  ```markdown
  The previous implementation used dynamic polymorphism with `std::unique_ptr` wrappers.
  To achieve optimal computational performance, we refactor the state into regular value types.
  ```

---

## 5. Release Notes & Process Standard

The single source of truth for the release lifecycle, pre-release checklist, and scientific release notes formatting is:
📄 **[docs/Release_Process.md](docs/Release_Process.md)**

When authoring release notes:
1. Use an objective, impersonal tone (no first-person pronouns or conversational greetings).
2. Follow the standard section layout defined in [docs/Release_Process.md](docs/Release_Process.md).
3. Ensure all sentence length limits (Rule 1) and British English spelling (Section 3) are strictly followed.

---

## 6. Diagramming Policy: Avoid Unnecessary Mermaid Diagrams

Do not use Mermaid diagrams for simple, linear, or text-first concepts.

### Principles:
1. **Prefer Standard Markdown:** Use bulleted lists, numbered steps, comparison tables, or ASCII / code blocks instead of Mermaid. Standard Markdown renders reliably across all browsers, mobile devices, diff viewers, and sandboxed environments without text clipping or layout bugs.
2. **Never Use Horizontal Flowcharts (`flowchart LR`):** Wide multi-stage horizontal diagrams collapse into unreadable narrow strips on GitHub's viewport.
3. **Avoid Text-Heavy Nodes:** Bounding-box calculation differences between system fonts and SVG containers cause persistent character truncation on GitHub.
4. **Reserve Diagrams for Genuine Need:** Only use diagrams when communicating non-trivial, multi-branch network topologies or complex state machines that cannot be understood as text. If a diagram is strictly required, use vertical orientation (`flowchart TD`) with minimal node labels.

---

## 7. Code Commenting & Mathematical Documentation Standards

RFL enforces a strict separation between public API documentation, internal implementation comments, and unit tests.

### 7.1 Separation by Codebase Layer

| Layer | File Types | Allowed Comment Format | Math Notation | Doxygen Tags Allowed? |
| :--- | :--- | :--- | :--- | :--- |
| **Public C++ Headers** | `src/core/**/*.hpp`, `include/**/*.hpp` | `/** ... */` (Javadoc style) | Clean ASCII / Markdown in `@brief`; LaTeX `\f$` in body | **Yes** (`@brief`, `@param`, `@return`, `@note`) |
| **Internal C++ Implementation** | `src/core/**/*.cpp` | `//` or `/* ... */` | Plain text, ASCII, Markdown backticks | **No** (Forbidden) |
| **Unit & Verification Tests** | `src/core/tests/**/*.cpp`, `tests/**/*.cpp` | `//` or `/* ... */` | Plain text, ASCII, Markdown (`$ ... $`) | **No** (Forbidden; breaks IDE hover cards) |
| **Python Bindings & Modules** | `python/rfl/`, `pyrfl` | `""" ... """` (Google / NumPyDoc) | Markdown code spans / LaTeX in Notes | **No** (Use Sphinx `:math:` only in Notes) |

### 7.2 The "WHY vs WHAT" Principle
Assume the reader understands C++ and Python syntax. Never restate syntax in English.
* ❌ **Do not write (WHAT):** `// Loops over all geometries and calls compareDelta.`
* ✅ **Write (WHY):** `// Verifies that local trace variations match brute-force recalculation across all admissible Euclidean and Lorentzian signatures.`
* ❌ **Do not write (WHAT):** `// Multiplies the diagonal entry by 2.`
* ✅ **Write (WHY):** `// The symmetric move δM = z·E_{ij} + z̄·E_{ji} collapses to z + z̄ = 2z on the diagonal (with z real).`

### 7.3 Dual-Tier Pattern for Public C++ Headers
Language servers (`clangd`, CLion) convert `@brief` into IDE hover cards. Raw LaTeX (`\f$`) renders as unreadable escape sequences.
* **Tier 1 (`@brief`):** Use plain ASCII and code spans so tooltips render cleanly in VS Code.
* **Tier 2 (Body & `@details`):** Use formal Doxygen LaTeX (`\f$ ... \f$` or `\f[ ... \f]`) for generated web manuals.

```cpp
/**
 * @class BarrettGlaserAction
 *
 * @brief Implements the Barrett-Glaser action: S(D) = g_2 Tr(D^2) + g_4 Tr(D^4).
 *
 * Evaluates the spectral action functional \f$ S(D) = g_2 \text{Tr}(D^2) + g_4 \text{Tr}(D^4) \f$.
 * Provides total action evaluation, local analytic trace variations, and gradients for HMC.
 *
 * Reference: Barrett & Glaser (2016), arXiv:1510.01377, Eq. (1.1).
 */
```

### 7.4 The 4-Part Invariant Framework for Unit Tests
Every test asserting mathematical, physical, or geometric invariants must document:
1. **Invariant Definition:** What symmetry or law is being tested (Hermiticity, gauge invariance, trace conservation)?
2. **Authoritative Literature Reference:** Author, year, paper title/arXiv ID, and equation or theorem number.
3. **Dual-Path Verification Mechanism:** Which optimized calculation is compared against what brute-force baseline?
4. **Tolerance Rationale:** Why this tolerance was chosen, based on machine epsilon ($\epsilon_{\text{mach}}$), matrix dimension ($M$), or cancellation error.

```cpp
// Applies an elementary variation delta_M to matrix x at entry (row, col).
//
// 1. Invariant:
//    Hermiticity preservation: (M + delta_M)† == M + delta_M.
//
// 2. Mathematical Move:
//    The Barrett-Glaser analytic trace variations evaluate the change under:
//    delta_M = z * E_{row,col} + conj(z) * E_{col,row}
//    - Off-diagonal (row != col): updates both (row, col) and (col, row).
//    - Diagonal (row == col): collapses to z + conj(z) = 2z (with Im(z) = 0),
//      scaling the diagonal update by 2.
//
// 3. Dual-Path Verification:
//    Used to perturb a copied Dirac operator to verify analytic delta formulas.
```

### 7.5 Tolerance Derivation Policy
Never assert bare magic constants without explaining their derivation:
* ❌ `EXPECT_NEAR(diff, delta, 1e-7);`
* ✅
  ```cpp
  // Scale-aware bound: accounts for matrix dimension M and double-well cancellation (EP-4).
  const double tol = 10.0 * M * std::numeric_limits<double>::epsilon() *
                     (std::abs(s_f) + std::abs(s_i) + trace_norm);
  EXPECT_NEAR(diff, delta, tol);
  ```

### 7.6 Automated Enforcement
The script `scripts/verify_comments.py` runs under `mise run lint`. It fails if:
* Doxygen tags (`\f$`, `\f[`, `/// \brief`, `/**`) appear in any test file (`*/tests/*`).
* Doxygen docblocks (`/**`) appear in implementation `.cpp` files.
* Raw LaTeX formulas appear in header `@brief` lines.

