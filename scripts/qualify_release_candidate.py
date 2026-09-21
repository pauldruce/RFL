#!/usr/bin/env python3
"""Executes pre-release qualification tasks for release candidates."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def run_cmd(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> tuple[int, str, str]:
    """Runs a shell command and returns returncode, stdout, and stderr."""
    res = subprocess.run(
        cmd,
        cwd=cwd or REPO_ROOT,
        capture_output=True,
        text=True,
        env=env or os.environ.copy(),
    )
    return res.returncode, res.stdout, res.stderr


def qualify_pypi(version: str, scratch_dir: Path) -> tuple[bool, str]:
    """Validates pyrfl pre-release package from PyPI in an isolated environment."""
    venv_dir = scratch_dir / "venv"
    if venv_dir.exists():
        shutil.rmtree(venv_dir)

    # 1. Create virtual environment
    code, out, err = run_cmd([sys.executable, "-m", "venv", str(venv_dir)])
    if code != 0:
        return False, f"Virtual environment creation failed: {err.strip()}"

    pip_bin = venv_dir / "bin" / "pip"
    python_bin = venv_dir / "bin" / "python"
    if not pip_bin.exists():
        pip_bin = venv_dir / "Scripts" / "pip.exe"
        python_bin = venv_dir / "Scripts" / "python.exe"

    # 2. Install pre-release package from PyPI
    clean_ver = version.lstrip("v")
    code, out, err = run_cmd([str(pip_bin), "install", "--quiet", "--pre", f"pyrfl=={clean_ver}"])
    if code != 0:
        return False, f"Failed to install pyrfl=={clean_ver} from PyPI: {err.strip()}"

    # 3. Execute smoke assertions in Python
    test_code = f"""
import rfl
assert rfl.__version__ == "{clean_ver}", f"Version mismatch: {{rfl.__version__}}"

dirac = rfl.DiracOperator(1, 3, 6)
assert dirac.get_matrix_dimension() > 0, "Dimension must be positive"

metro = rfl.Metropolis(-1.0, 1.0, 0.05, 50, 42)
acceptance = metro.update_dirac(dirac)
assert 0.0 <= acceptance <= 1.0, f"Acceptance rate out of range: {{acceptance}}"

eigenvalues = dirac.get_eigenvalues()
assert len(eigenvalues) > 0, "Eigenvalue spectrum must not be empty"
print(f"Verified pyrfl=={{rfl.__version__}}, eigenvalues={{len(eigenvalues)}}, acceptance={{acceptance*100:.2f}}%")
"""
    code, out, err = run_cmd([str(python_bin), "-c", test_code])
    if code != 0:
        return False, f"Python smoke verification failed: {err.strip()}"

    # 4. Execute examples/python/main.py with installed package
    python_example = REPO_ROOT / "examples" / "python" / "main.py"
    if python_example.exists():
        code, out, err = run_cmd([str(python_bin), str(python_example)])
        if code != 0:
            return False, f"examples/python/main.py failed with installed package: {err.strip()}"

    return True, f"pyrfl=={clean_ver} installed from PyPI and validated successfully."


def qualify_fetchcontent(tag: str, scratch_dir: Path) -> tuple[bool, str]:
    """Validates downstream CMake consumption and target isolation via FetchContent."""
    canary_source_dir = REPO_ROOT / "examples" / "downstream_canary"
    consumer_dir = scratch_dir / "consumer"
    if consumer_dir.exists():
        shutil.rmtree(consumer_dir)

    # Use canonical in-tree downstream canary project
    shutil.copytree(canary_source_dir, consumer_dir)

    # 1. Configure CMake with specified Git tag
    build_dir = consumer_dir / "build"
    code, out, err = run_cmd(
        ["cmake", "-B", str(build_dir), "-S", str(consumer_dir), f"-DRFL_GIT_TAG={tag}"],
        cwd=consumer_dir,
    )
    if code != 0:
        return False, f"Consumer CMake configure failed: {err.strip()}"

    # 2. Verify target isolation: internal RFL test targets must NOT be present
    ninja_file = build_dir / "build.ninja"
    makefile = build_dir / "Makefile"
    build_spec = ""
    if ninja_file.exists():
        build_spec = ninja_file.read_text(encoding="utf-8", errors="ignore")
    elif makefile.exists():
        build_spec = makefile.read_text(encoding="utf-8", errors="ignore")

    if "rfl_tests" in build_spec or "rfl_performance_tests" in build_spec:
        return False, "Target isolation failure: internal RFL test targets leaked into consumer build."

    # 3. Build consumer application
    code, out, err = run_cmd(
        ["cmake", "--build", str(build_dir), "--target", "downstream_canary"],
        cwd=consumer_dir,
    )
    if code != 0:
        return False, f"Consumer CMake build failed: {err.strip()}"

    # 4. Run consumer application
    binary_path = build_dir / "downstream_canary"
    if not binary_path.exists():
        binary_path = build_dir / "downstream_canary.exe"

    code, out, err = run_cmd([str(binary_path)], cwd=consumer_dir)
    if code != 0:
        return False, f"Consumer executable failed to run: {err.strip()}"

    return True, f"Downstream FetchContent ({tag}) built and verified with target isolation."


def qualify_examples() -> tuple[bool, str]:
    """Validates local CMake, Makefile, and direct compiler example workflows."""
    # 1. CMake example targets
    example_targets = ["main", "downstream_canary"]
    candidates = [
        REPO_ROOT / "build" / "examples" / "cpp" / "main",
        REPO_ROOT / "build" / "examples" / "downstream_canary" / "downstream_canary",
    ]

    # Detect whether optional legacy targets exist (requires GSL)
    code_check, help_out, _ = run_cmd(["cmake", "--build", "build", "--target", "help"], cwd=REPO_ROOT)
    if code_check == 0 and "mauro_thesis_mmc" in help_out:
        example_targets.extend(["mauro_thesis_mmc", "hmc_tuning"])
        candidates.extend([
            REPO_ROOT / "build" / "examples" / "case_studies" / "mauro_thesis_mmc" / "mauro_thesis_mmc",
            REPO_ROOT / "build" / "examples" / "case_studies" / "hmc_tuning" / "hmc_tuning",
        ])

    code, out, err = run_cmd(
        ["cmake", "--build", "build", "--target"] + example_targets,
        cwd=REPO_ROOT,
    )
    if code != 0:
        return False, f"CMake example targets failed to build: {err.strip()}"

    for candidate_path in candidates:
        target_bin = candidate_path
        if not target_bin.exists():
            target_bin = target_bin.with_suffix(".exe")
        if not target_bin.exists():
            flat_candidate = candidate_path.parent.parent / candidate_path.name
            if flat_candidate.exists() or flat_candidate.with_suffix(".exe").exists():
                target_bin = flat_candidate if flat_candidate.exists() else flat_candidate.with_suffix(".exe")
        # Run within target directory so generated files (e.g. example_S.txt) stay within build/
        code, out, err = run_cmd([str(target_bin)], cwd=target_bin.parent)
        if code != 0:
            return False, f"Example executable {target_bin.name} failed: {err.strip()}"

    # 2. Standalone Makefile
    makefile_dir = REPO_ROOT / "examples" / "cpp"
    run_cmd(["make", "clean"], cwd=makefile_dir)
    try:
        code, out, err = run_cmd(["make", "run"], cwd=makefile_dir)
        if code != 0:
            return False, f"examples/cpp Makefile workflow failed: {err.strip()}"
    finally:
        run_cmd(["make", "clean"], cwd=makefile_dir)

    # 3. Direct GCC / Clang script
    script_path = REPO_ROOT / "examples" / "cpp" / "compile_gcc.sh"
    if script_path.exists():
        code, out, err = run_cmd([str(script_path)], cwd=REPO_ROOT)
        if code != 0:
            return False, f"compile_gcc.sh workflow failed: {err.strip()}"

    return True, "All example workflows (CMake, Makefile, direct compiler) passed."


def get_latest_git_tag() -> str | None:
    """Returns the most recent Git release tag."""
    try:
        res = subprocess.run(
            ["git", "tag", "--sort=-v:refname", "--list", "v[0-9]*"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        tags = [t.strip() for t in res.stdout.splitlines() if t.strip()]
        return tags[0] if tags else None
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Pre-Release Qualification Test Suite for RFL")
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Release candidate tag to qualify (e.g. v0.2.0rc2, defaults to latest tag)",
    )
    parser.add_argument(
        "--skip-pypi",
        action="store_true",
        default=False,
        help="Skip PyPI package installation and validation check",
    )
    parser.add_argument(
        "--skip-fetchcontent",
        action="store_true",
        default=False,
        help="Skip downstream FetchContent CMake validation check",
    )
    parser.add_argument(
        "--skip-examples",
        action="store_true",
        default=False,
        help="Skip example application workflow checks",
    )
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        default=None,
        help="Custom scratch directory for qualification artifacts (defaults to system temporary directory)",
    )
    parser.add_argument(
        "--keep-artifacts",
        action="store_true",
        default=False,
        help="Retain qualification artifacts in build/qualification_<tag> instead of deleting temporary files",
    )

    args = parser.parse_args()
    tag = args.tag or get_latest_git_tag()
    if not tag:
        print("❌ Error: No release tag specified or found in Git history.", file=sys.stderr)
        return 1

    print("============================================================")
    print(f"🔎 RFL Pre-Release Qualification Suite: {tag}")
    print("============================================================\n")

    temp_dir_ctx: tempfile.TemporaryDirectory | None = None
    if args.scratch_dir:
        scratch_base = args.scratch_dir.resolve()
        scratch_base.mkdir(parents=True, exist_ok=True)
    elif args.keep_artifacts:
        scratch_base = REPO_ROOT / "build" / f"qualification_{tag.replace('.', '_')}"
        scratch_base.mkdir(parents=True, exist_ok=True)
    else:
        temp_dir_ctx = tempfile.TemporaryDirectory(prefix=f"rfl_qualification_{tag.replace('.', '_')}_")
        scratch_base = Path(temp_dir_ctx.name)

    results: list[tuple[str, bool, str]] = []

    try:
        # Task 1: PyPI Verification
        if not args.skip_pypi:
            print(f"[1/3] Qualifying PyPI package (pyrfl=={tag.lstrip('v')})...")
            passed, msg = qualify_pypi(tag, scratch_base)
            results.append(("PyPI Pre-Release Wheel", passed, msg))
            print(f"      {'✅' if passed else '❌'} {msg}\n")
        else:
            print("[1/3] Skipping PyPI qualification check.\n")

        # Task 2: FetchContent Consumer Verification
        if not args.skip_fetchcontent:
            print(f"[2/3] Qualifying downstream CMake FetchContent ({tag})...")
            passed, msg = qualify_fetchcontent(tag, scratch_base)
            results.append(("Downstream FetchContent", passed, msg))
            print(f"      {'✅' if passed else '❌'} {msg}\n")
        else:
            print("[2/3] Skipping FetchContent qualification check.\n")

        # Task 3: Local Examples Verification
        if not args.skip_examples:
            print("[3/3] Qualifying example workflows (CMake, Makefile, direct)...")
            passed, msg = qualify_examples()
            results.append(("Example Application Workflows", passed, msg))
            print(f"      {'✅' if passed else '❌'} {msg}\n")
        else:
            print("[3/3] Skipping example application qualification checks.\n")
    finally:
        if temp_dir_ctx is not None:
            temp_dir_ctx.cleanup()

    # Summary Report
    print("============================================================")
    print(f"📊 Qualification Summary Report: {tag}")
    print("============================================================")
    all_passed = True
    for name, passed, detail in results:
        status_str = "PASS" if passed else "FAIL"
        symbol = "✅" if passed else "❌"
        print(f"{symbol} [{status_str:4s}] {name:30s} : {detail}")
        if not passed:
            all_passed = False

    print("============================================================")
    if all_passed:
        print(f"🎉 Candidate '{tag}' PASSED all pre-release qualification tasks.")
        print("   Ready for promotion to final release.")
        return 0
    else:
        print(f"❌ Candidate '{tag}' FAILED qualification tasks.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
