# RFL Release Process & Scientific Release Notes Guide

This document establishes the official release lifecycle, pre-release checklist, and release notes writing standard for the `RFL` (Random Fuzzy Library) project.

---

## 1. Release Philosophy & Governance

1. **Semantic Versioning & Beta Lifecycle:** Releases follow `vMAJOR.MINOR.PATCH` (e.g. `v0.5.0`).
   - **Beta Development Phase (`v0.y.z`):** The library is currently in an active research and architectural modernisation phase. Under Semantic Versioning rules, minor version increments (`v0.5.0` → `v0.6.0`) may introduce breaking API changes or refactors while the architecture evolves toward `v1.0.0`.
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
   - Transitioning from draft to published initiates the automated wheel build and PyPI distribution pipeline.

---

## 2. Pre-Releases & Release Candidates (RCs)

Following the convention of major scientific libraries (such as NumPy, SciPy, PyTorch, and Armadillo), significant releases use **Release Candidates** (e.g. `v0.5.0rc1`):

```mermaid
flowchart TD
    subgraph "Phase 1: Pre-Release (RC)"
        A["Tag: v0.5.0rc1\n(git tag v0.5.0rc1)"] --> B["Create Draft Pre-Release\n(gh release create --draft --prerelease)"]
        B --> C["Review Draft Notes in UI\n(Zero workflows triggered)"]
        C --> D["Publish Pre-Release\n(gh release edit --draft=false)"]
        D --> E["Automated CI/CD\n- Build wheels & sdist\n- Upload assets\n- Publish to PyPI as pre-release"]
        E --> F["Testing Period (24–72h)\n- pip install --pre rfl\n- CMake FetchContent (GIT_TAG v0.5.0rc1)"]
    end

    subgraph "Phase 2: Final Promotion"
        F --> G{"Validation successful?"}
        G -- Issues found --> H["Fix on branch → Tag v0.5.0rc2"]
        H --> B
        G -- Clean --> I["Tag Final: v0.5.0\n(git tag v0.5.0)"]
        I --> J["Create Draft Release\n(gh release create --draft)"]
        J --> K["Review Draft Notes in UI\n(Zero workflows triggered)"]
        K --> L["Publish Final Release\n(Official PyPI default & Latest tag)"]
    end
```

### 2.1 How Pre-Releases Work Across Ecosystems

1. **Naming Standard (PEP 440 & Git SemVer):**
   - Use `vX.Y.Zrc1` (e.g. `v0.5.0rc1`).
   - This tag format is natively recognised by Git, GitHub, `pip`, and `scikit-build-core`.
2. **GitHub Releases Behaviour:**
   - Pre-releases are flagged with `--prerelease` (or the "Set as a pre-release" checkbox).
   - GitHub displays a `Pre-release` badge and retains the previous release as `Latest`.
3. **PyPI & `pip` Behaviour:**
   - PyPI automatically marks `0.5.0rc1` as a pre-release.
   - A standard `pip install rfl` will **never** install a pre-release by default.
   - Downstream researchers must explicitly opt in with:
     ```bash
     pip install --pre rfl
     # or
     pip install rfl==0.5.0rc1
     ```
4. **C++ & CMake `FetchContent` Behaviour:**
   - Downstream C++ solvers test the release candidate by pinning the RC git tag:
     ```cmake
     FetchContent_Declare(
         RFL
         GIT_REPOSITORY https://github.com/pauldruce/RFL.git
         GIT_TAG        v0.5.0rc1
     )
     ```

---

## 3. The Complete Release Lifecycle

### Step 1: Pre-Release Checklist
Before tagging any release or candidate:
* Ensure all CI workflows on `main` are green.
* Verify local builds and tests pass (`ctest`, `pytest`).

### Step 2: Milestone Triage & EP Status
* Verify that all GitHub Issues and PRs for the milestone are merged and closed.
* Update the relevant Enhancement Proposal in [docs/eps/](eps/) to reflect current milestone status.

### Step 3: Tag, Draft, and Publish a Release Candidate (`vX.Y.Zrc1`)
1. Create and push the release candidate tag:
   ```bash
   git checkout main
   git pull origin main
   git tag v0.5.0rc1
   git push origin v0.5.0rc1
   ```
2. Create the Draft Pre-Release using the committed release notes:
   ```bash
   gh release create v0.5.0rc1 --draft --prerelease --title "v0.5.0rc1: Release Candidate 1" --notes-file docs/releases/v0.5.0.md
   ```
3. Review the rendered release notes in the GitHub Web UI.
4. Publish the pre-release:
   ```bash
   gh release edit v0.5.0rc1 --draft=false
   ```
   The release workflow compiles wheels, packages `sdist`, attaches release assets, and uploads the pre-release to PyPI.

### Step 4: Release Candidate Testing
During the testing window (24–72 hours for `0.y.z` releases):
1. **Python Wheels:** Verify in a clean virtual environment:
   ```bash
   python -m venv test_env && source test_env/bin/activate
   pip install --pre rfl
   python -c "import rfl; print(rfl.__version__)"
   ```
2. **C++ CMake Integration:** Verify `FetchContent` in an external test project linking `RFL::core`.

### Step 5: Tag, Draft, and Publish Final Release (`vX.Y.Z`)
When the release candidate is validated with zero critical defects:
1. Tag the release commit:
   ```bash
   git tag v0.5.0
   git push origin v0.5.0
   ```
2. Create the Draft Release:
   ```bash
   gh release create v0.5.0 --draft --title "v0.5.0: Multi-Platform Package and Binary Distribution" --notes-file docs/releases/v0.5.0.md
   ```
3. Review the rendered release notes in the GitHub Web UI.
4. Publish the final release:
   ```bash
   gh release edit v0.5.0 --draft=false
   ```
   The final packages become the default on PyPI (`pip install rfl`), assets are attached to GitHub Releases, and the milestone is closed.

---

## 4. In-Tree Release Notes Storage (`docs/releases/`)

Following scientific software conventions (NumPy, SciPy, Boost), all release notes are committed directly to the repository under **`docs/releases/vX.Y.Z.md`**.

### 4.1 Canonical Template
The single source of truth for authoring release notes is:
📄 **[docs/release-note-template.md](release-note-template.md)**

To author release notes for a new version:
1. Copy `docs/release-note-template.md` to `docs/releases/vX.Y.Z.md`.
2. Populate the breaking changes alert box, highlights, and component changes.
3. Commit `docs/releases/vX.Y.Z.md` on the release PR for peer review.
4. Pass `--notes-file docs/releases/vX.Y.Z.md` when creating the GitHub Release.

### 4.2 Historical Pre-Releases (`v0.1.0` – `v0.4.0`)
Historical development milestones are backfilled in `docs/releases/` (`v0.1.0.md` through `v0.4.0.md`). These milestones were unpackaged internal iterations. For all Python and modern C++ research workflows, use **RFL `v0.5.0` or later**.

---

## 5. Testing the Release Pipeline (Dry Run)

To test wheel compilation and packaging without uploading to PyPI or modifying release assets, trigger a manual dry run:

```bash
gh workflow run release.yml -f dry_run=true -f publish_to_pypi=false
```
