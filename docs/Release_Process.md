# RFL Release Process & Scientific Release Notes Guide

This document establishes the official release lifecycle, pre-release checklist, and release notes writing standard for the `RFL` (Random Fuzzy Library) project.

---

## 1. Release Philosophy & Governance

1. **Semantic Versioning & Beta Lifecycle:** Releases follow `vMAJOR.MINOR.PATCH` (e.g. `v0.2.0`).
   - **Beta Development Phase (`v0.y.z`):** The library is currently in an active research phase. Minor version increments (`v0.1.0` → `v0.2.0`) may introduce breaking API changes before `v1.0.0`.
   - **Stable Production Phase (`v1.0.0+`):** After `v1.0.0`, breaking changes occur only across MAJOR version increments, with a formal deprecation period across MINOR releases.
2. **Controlled Language (ASD-STE100 & British English):**
   - Release documentation must follow controlled vocabulary defined in [docs/Glossary.md](Glossary.md).
   - Use British English spelling (*standardisation*, *optimisation*, *behaviour*, *modelling*).
3. **Impersonal Scientific Tone:**
   - Write in an objective, third-person perspective.
   - Do not use first-person pronouns (*"we"*, *"I"*, *"our"*).
   - Avoid conversational openings (e.g. *"We are pleased to announce..."*).
4. **Draft-First Review Gate:**
   - Always create releases and release candidates as **Drafts** first (`--draft`).
   - Review rendered release notes, breaking change warnings, and links in the GitHub UI before publishing.
   - Publishing initiates the automated wheel build and PyPI distribution pipeline.

---

## 2. Pre-Releases & Release Candidates (RCs)

Following major scientific library conventions (such as NumPy and SciPy), significant releases use **Release Candidates** (e.g. `v0.2.0rc1`):

| Stage | Action | Command / Trigger | Automation & Verification |
| :--- | :--- | :--- | :--- |
| **1. Tag RC** | Create release candidate tag | `git tag v0.2.0rc1 && git push origin v0.2.0rc1` | Records candidate commit point |
| **2. Draft Pre-Release** | Author notes for review | `gh release create v0.2.0rc1 --draft --prerelease` | Zero CI workflows triggered while in draft |
| **3. Publish RC** | Trigger packaging pipeline | `gh release edit v0.2.0rc1 --draft=false` | Builds wheels, uploads assets, publishes pre-release to PyPI |
| **4. Test & Qualify** | Testing period (24–72h) | `pip install --pre pyrfl` | Downstream verification on user machines |
| **5. Resolve or Promote** | Fix defects or promote to final | `git tag v0.2.0 && gh release create v0.2.0` | Official release on PyPI and GitHub `Latest` tag |

### 2.1 How Pre-Releases Work Across Ecosystems

1. **Naming Standard (PEP 440 & Git SemVer):**
   - Use `vX.Y.Zrc1` (e.g. `v0.2.0rc1`).
   - This tag format is natively recognised by Git, GitHub, `pip`, and `scikit-build-core`.
2. **GitHub Releases Behaviour:**
   - Pre-releases are flagged with `--prerelease` (or the "Set as a pre-release" checkbox).
   - GitHub displays a `Pre-release` badge and retains the previous release as `Latest`.
3. **PyPI & `pip` Behaviour:**
   - PyPI automatically marks `0.2.0rc1` as a pre-release.
   - Standard `pip install pyrfl` will never install a pre-release by default.
   - Downstream researchers must explicitly opt in:
     ```bash
     pip install --pre pyrfl
     # or
     pip install pyrfl==0.2.0rc1
     ```
4. **C++ & CMake `FetchContent` Behaviour:**
   - Downstream C++ solvers test the release candidate by pinning the RC git tag:
     ```cmake
     FetchContent_Declare(
         RFL
         GIT_REPOSITORY https://github.com/pauldruce/RFL.git
         GIT_TAG        v0.2.0rc1
     )
     ```

---

## 3. The Complete Release Lifecycle

### Step 1: Pre-Release Checklist
Before tagging any release or candidate:
* Ensure all CI workflows on `main` pass.
* Verify local builds and tests pass (`ctest`, `pytest`).

### Step 2: Milestone Triage & EP Status
* Verify that all GitHub issues and PRs for the milestone are merged and closed.
* Update the relevant Enhancement Proposal in [docs/eps/](eps/) to reflect current milestone status.

### Step 3: Tag, Draft, and Publish a Release Candidate (`vX.Y.Zrc1`)
1. Update `CHANGELOG.md` for the target version and verify with:
   ```bash
   python3 scripts/verify_changelog.py --tag vX.Y.Zrc1
   ```
2. Create and push the release candidate tag:
   ```bash
   git checkout main
   git pull origin main
   git tag vX.Y.Zrc1
   git push origin vX.Y.Zrc1
   ```
3. Create the Draft Pre-Release on GitHub:
   ```bash
   gh release create vX.Y.Zrc1 --draft --prerelease --title "vX.Y.Zrc1: Release Candidate 1" --generate-notes
   ```
4. Review the rendered release notes in the GitHub Web UI.
5. Publish the pre-release:
   ```bash
   gh release edit vX.Y.Zrc1 --draft=false
   ```
   Publishing triggers wheel compilation, sdist packaging, and PyPI pre-release deployment.

### Step 4: Pre-Release Qualification Protocol
During the qualification window (24–72 hours for `0.y.z` releases), execute the automated qualification suite:

```bash
python3 scripts/qualify_release_candidate.py --tag vX.Y.ZrcN
```

The qualification suite executes three rigorous integration audits:

1. **Task 1: PyPI Pre-Release Package Verification (`pyrfl`):**
   * Installs the pre-release package from PyPI into a clean virtual environment (`pip install --pre pyrfl==X.Y.ZrcN`).
   * Validates runtime version strings, Dirac operator initialisation, Metropolis Monte Carlo sweeps, and eigenvalue extraction.
   * Executes the user-facing example application [`examples/python/main.py`](../examples/python/main.py) against the installed package.

2. **Task 2: Downstream C++ Consumer Verification (`FetchContent`):**
   * Configures an isolated downstream CMake project that consumes `RFL::core` via `FetchContent`.
   * Verifies downstream target isolation: internal RFL test targets (`rfl_tests`, `rfl_performance_tests`, `googletest`) must not leak into consumer targets.
   * Compiles and executes the consumer binary successfully.

3. **Task 3: Local C++ Example Workflow Verification:**
   * Compiles and executes all CMake example binaries (`main`, `mauro_thesis_mmc`, `hmc_tuning`).
   * Tests the standalone Makefile workflow in `examples/cpp` (`make`, `make run`, `make clean`).
   * Tests direct compiler compilation via [`examples/cpp/compile_gcc.sh`](../examples/cpp/compile_gcc.sh).

### Step 5: Tag, Draft, and Publish Final Release (`vX.Y.Z`)
When the release candidate completes verification without critical defects:
1. Verify `CHANGELOG.md` with:
   ```bash
   python3 scripts/verify_changelog.py --tag vX.Y.Z
   ```
2. Tag the release commit:
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
3. Create the Draft Release on GitHub:
   ```bash
   gh release create vX.Y.Z --draft --title "vX.Y.Z: Release Title" --generate-notes
   ```
4. Review rendered release notes in the GitHub Web UI.
5. Publish the final release:
   ```bash
   gh release edit vX.Y.Z --draft=false
   ```
   Publishing promotes packages to the PyPI default, attaches release assets, and closes the milestone.

---

## 4. Release History & Changelog Governance (`CHANGELOG.md`)

RFL maintains release notes and version history directly in **`CHANGELOG.md`** following the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) specification.

### 4.1 Automated Changelog Verification Gate
The automated test script `scripts/verify_changelog.py` validates that:
* Every release tag has a corresponding version section in `CHANGELOG.md`.
* All pull requests merged since the previous release are documented.
* No unresolved placeholder tokens (`TODO`, `TBD`) exist.

Run the verification gate locally:
```bash
python3 scripts/verify_changelog.py --tag vX.Y.Z
```

### 4.2 Initial Release Baseline (`v0.1.0`)
RFL versioning formally begins with `v0.1.0` as the initial packaged release with binary wheels and CMake target exports.

### 4.3 Patch Releases & Maintenance Branch Workflow
When a bug fix patch (`vX.Y.1`) is needed while `main` develops future versions (`vX.(Y+1).0`), use the Maintenance Branch Workflow:
1. **Fix on `main` First:** Land bug fixes on `main` through a pull request to prevent regressions.
2. **Backport to Maintenance Branch:** Cherry-pick the bug fix commit to `maintenance/X.Y.x`.
3. **Update Changelog:** Update `CHANGELOG.md` on the maintenance branch to document the repaired issues. Patch releases must never contain breaking changes.
4. **Tag & Publish:** Tag `vX.Y.1` from the maintenance branch and publish using `gh release create`.
5. **Forward-Port Notes to `main`:** Merge changelog updates back to `main` to retain full release history.

---

## 5. Testing the Release Pipeline & TestPyPI Qualification

### 5.1 Automated Pull Request Testing
Any pull request that modifies `.github/workflows/release.yml`, `pyproject.toml`, or `src/python_bindings/**` automatically executes the complete packaging matrix. The workflow compiles binary wheels across Linux, macOS, and Windows, and executes the `pytest` test suite without publishing assets.

### 5.2 TestPyPI Publication & Pre-Flight Qualification
To publish packages to TestPyPI (`test.pypi.org`) for downstream verification before an official release, dispatch the release workflow with `publish_to_testpypi=true`:

```bash
gh workflow run release.yml -f publish_to_testpypi=true -f publish_to_pypi=false
```

Once publication completes, verify installation in a clean environment:

```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ --pre pyrfl
python -c "import rfl; print('RFL loaded successfully from TestPyPI:', rfl)"
```

### 5.3 Dry Run (No Uploads)
To test wheel compilation and packaging without publishing to any index or repository, trigger a manual dry run:

```bash
gh workflow run release.yml -f dry_run=true -f publish_to_pypi=false -f publish_to_testpypi=false
```
