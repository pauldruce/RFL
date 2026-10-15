#!/usr/bin/env python3
"""scripts/verify_tdd_commits.py

Atomic Two-Commit Red-Green Protocol Inspector.
Verifies test-driven development compliance for feature pull requests under EP-5 (REQ-TRC-007).

Resolves the Red Anchor deterministically as the immediate parent of the first feat: commit.
Checks that stubs fail on assertions (<failure>) at Red, and pass (<status="PASSED">) at Green.
Uses in-place checkout on Ubuntu runners to preserve ccache and Ninja build state.

Conforms to EP-5 (REQ-TRC-007).
Uses standard Python library only (zero external dependencies).
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


FEAT_COMMIT_PATTERN = re.compile(r"^feat(\([^)]+\))?:", re.IGNORECASE)


def resolve_red_and_green_anchors(
    base_ref: str, repo_root: Path
) -> tuple[str | None, str | None, list[str]]:
    """Resolve the Red Anchor SHA and Green Anchor SHA from Git history.

    The Red Anchor is the immediate parent of the first feat: commit.
    This captures the Red state after all test stubs are registered.
    The Green Anchor is the current HEAD.
    """
    # Query commits from base_ref to HEAD in chronological order with parents
    cmd = [
        "git",
        "rev-list",
        "--reverse",
        "--parents",
        f"{base_ref}..HEAD",
    ]
    try:
        process_result = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as err:
        return None, None, [f"Failed to query git rev-list for {base_ref}..HEAD: {err}"]

    lines = [line.strip() for line in process_result.stdout.splitlines() if line.strip()]
    if not lines:
        return None, None, ["No commits found between base and HEAD."]

    first_feat_sha: str | None = None
    red_anchor_sha: str | None = None

    for line in lines:
        parts = line.split()
        commit_sha = parts[0]
        parent_shas = parts[1:]

        # Query commit subject message
        subject_cmd = ["git", "log", "-1", "--format=%s", commit_sha]
        subject_result = subprocess.run(subject_cmd, cwd=repo_root, capture_output=True, text=True, check=True)
        subject = subject_result.stdout.strip()

        if FEAT_COMMIT_PATTERN.match(subject):
            first_feat_sha = commit_sha
            if not parent_shas:
                return None, None, [f"First feature commit '{commit_sha}' has no parent commit."]
            red_anchor_sha = parent_shas[0]
            break

    if not first_feat_sha or not red_anchor_sha:
        return None, None, []

    # Green anchor is always HEAD
    head_result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo_root, capture_output=True, text=True, check=True)
    green_anchor_sha = head_result.stdout.strip()

    return red_anchor_sha, green_anchor_sha, []


def verify_tdd_protocol(
    base_ref: str,
    repo_root: Path,
    build_dir: Path | None = None,
    dry_run: bool = False,
    xml_output_dir: Path | None = None,
) -> int:
    """Execute the atomic Red-Green qualification protocol."""
    print("=" * 70)
    print(" RFL Atomic Red-Green Protocol Verification (EP-5 / REQ-TRC-007)")
    print("=" * 70)

    red_sha, green_sha, errors = resolve_red_and_green_anchors(base_ref, repo_root)
    if errors:
        print("\n❌ Commit topology error:")
        for err in errors:
            print(f"  • {err}")
        print("=" * 70)
        return 1

    if not red_sha or not green_sha:
        print("ℹ️  No 'feat:' commits found between base and HEAD.")
        print("   Treated as maintenance or refactoring PR (non-blocking pass).")
        print("=" * 70)
        return 0

    print(f"🔴 Red Anchor Resolved:   {red_sha}")
    print(f"🟢 Green Anchor Resolved: {green_sha}")

    if dry_run:
        print("\n✅ Red-Green commit topology successfully verified (dry run).")
        print("=" * 70)
        return 0

    # Ensure build directory exists for in-place checkout
    if build_dir is None:
        build_dir = repo_root / "build"

    if not build_dir.exists():
        print(f"❌ Build directory does not exist for in-place checkout: {build_dir}")
        print("   Run CMake configuration before executing TDD protocol verification.")
        return 1

    out_dir = xml_output_dir or build_dir
    red_xml = out_dir / "junit_red.xml"
    green_xml = out_dir / "junit_green.xml"

    # Save original HEAD ref to guarantee clean restore
    head_ref_result = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    original_ref = head_ref_result.stdout.strip()
    if original_ref == "HEAD":
        original_ref = green_sha

    step_summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    summary_lines: list[str] = [
        "### 🧪 Atomic Red-Green TDD Protocol Verification",
        "",
        f"- **Red Anchor Commit:** `{red_sha[:8]}` (Parent of first `feat:` commit)",
        f"- **Green Anchor Commit:** `{green_sha[:8]}` (`HEAD`)",
        "",
        "| Stage | Target Anchor | Expected Result | Actual Result | Verification |",
        "| :--- | :--- | :--- | :--- | :--- |",
    ]

    try:
        # 1. In-place checkout of Red Anchor
        print(f"\n🔄 Checking out Red Anchor: {red_sha}...")
        subprocess.run(["git", "checkout", red_sha], cwd=repo_root, check=True)

        print("🔨 Compiling test stubs at Red Anchor...")
        subprocess.run(["cmake", "--build", str(build_dir)], cwd=repo_root, check=True)

        print("🧪 Executing tests at Red Anchor (expecting assertion <failure>)...")
        ctest_red = subprocess.run(
            ["ctest", "--test-dir", str(build_dir), "--output-junit", str(red_xml)],
            cwd=repo_root,
            capture_output=True,
            text=True,
        )

        # Run verify_traceability with --red-anchor
        red_verify_cmd = [
            sys.executable,
            str(repo_root / "scripts" / "verify_traceability.py"),
            "--check-pr",
            "--base",
            base_ref,
            "--red-anchor",
            "--verify-xml",
            str(red_xml),
        ]
        red_result = subprocess.run(red_verify_cmd, cwd=repo_root, capture_output=True, text=True)

        if red_result.returncode != 0:
            print("\n❌ Red Anchor validation failed:")
            print(red_result.stdout)
            print(red_result.stderr)
            summary_lines.append(
                f"| **Red Phase** | `{red_sha[:8]}` | `<failure>` on stubs | Failed or uncaught error | ❌ FAILED |"
            )
            return 1

        print("✅ Red Anchor successfully verified (failed assertion on stubs, zero uncaught errors).")
        summary_lines.append(
            f"| **Red Phase** | `{red_sha[:8]}` | `<failure>` on stubs | Assertion `<failure>` confirmed | ✅ PASSED |"
        )

        # 2. In-place checkout of Green Anchor
        print(f"\n🔄 Checking out Green Anchor: {green_sha}...")
        subprocess.run(["git", "checkout", green_sha], cwd=repo_root, check=True)

        print("🔨 Incremental compilation at Green Anchor...")
        subprocess.run(["cmake", "--build", str(build_dir)], cwd=repo_root, check=True)

        print("🧪 Executing tests at Green Anchor (expecting all tests pass)...")
        subprocess.run(
            ["ctest", "--test-dir", str(build_dir), "--output-junit", str(green_xml)],
            cwd=repo_root,
            check=True,
        )

        green_verify_cmd = [
            sys.executable,
            str(repo_root / "scripts" / "verify_traceability.py"),
            "--check-pr",
            "--base",
            base_ref,
            "--verify-xml",
            str(green_xml),
        ]
        green_result = subprocess.run(green_verify_cmd, cwd=repo_root, capture_output=True, text=True)

        if green_result.returncode != 0:
            print("\n❌ Green Anchor validation failed:")
            print(green_result.stdout)
            print(green_result.stderr)
            summary_lines.append(
                f"| **Green Phase** | `{green_sha[:8]}` | All tests pass | Tests failing | ❌ FAILED |"
            )
            return 1

        print("✅ Green Anchor successfully verified (all tests passed).")
        summary_lines.append(
            f"| **Green Phase** | `{green_sha[:8]}` | All tests pass | 100% tests passed | ✅ PASSED |"
        )

    finally:
        # Always restore original checkout ref
        print(f"\n🔄 Restoring original working ref: {original_ref}...")
        subprocess.run(["git", "checkout", original_ref], cwd=repo_root, check=False)

        if step_summary_path:
            try:
                with open(step_summary_path, "a", encoding="utf-8") as f:
                    f.write("\n".join(summary_lines))
                    f.write("\n\n")
            except Exception:
                pass

    print("\n" + "=" * 70)
    print("✅ Atomic Red-Green Protocol successfully qualified.")
    print("=" * 70)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RFL Two-Commit Red-Green Protocol Inspector (EP-5 / REQ-TRC-007)"
    )
    parser.add_argument(
        "--base",
        type=str,
        default="origin/main",
        help="Git base reference (default: origin/main)",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Repository root directory",
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=None,
        help="Build directory for in-place checkout and incremental compilation",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Verify Red/Green commit topology without checking out or compiling",
    )
    parser.add_argument(
        "--xml-output-dir",
        type=Path,
        default=None,
        help="Directory to write JUnit XML test reports",
    )

    args = parser.parse_args()
    return verify_tdd_protocol(
        base_ref=args.base,
        repo_root=args.repo_root.resolve(),
        build_dir=args.build_dir.resolve() if args.build_dir else None,
        dry_run=args.dry_run,
        xml_output_dir=args.xml_output_dir.resolve() if args.xml_output_dir else None,
    )


if __name__ == "__main__":
    sys.exit(main())
