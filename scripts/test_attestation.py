#!/usr/bin/env python3
"""
test_attestation.py

Manages cryptographic tree-hash test attestation receipts for fast pre-push skips.
Follows ASD-STE100 principles with deterministic validation and zero external dependencies.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

MAX_ATTESTATION_AGE_SECONDS = 7200  # 2 hours


def get_git_dir() -> Path | None:
    """Return the absolute path to the .git directory (handles worktrees)."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--git-dir"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        return Path(out).resolve()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_current_state() -> tuple[str, str, bool] | None:
    """Return (head_sha, tree_sha, is_clean)."""
    try:
        head_sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()

        tree_sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD^{tree}"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()

        # Check if working directory or staged index has modifications
        diff_status = subprocess.run(
            ["git", "diff-index", "--quiet", "HEAD", "--"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        is_clean = (diff_status == 0)

        return head_sha, tree_sha, is_clean
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_attestation_file(target_dir: Path | None = None) -> Path:
    """Return path to test_attestation.json in target_dir, .private/, build/, or .git/."""
    if target_dir is not None:
        return target_dir / "test_attestation.json"
    repo_root = Path(__file__).resolve().parent.parent
    for d in [repo_root / ".private", repo_root / "build"]:
        if d.exists() and os.access(d, os.W_OK):
            return d / "test_attestation.json"
    git_dir = get_git_dir()
    if git_dir and os.access(git_dir, os.W_OK):
        return git_dir / "test_attestation.json"
    return repo_root / ".private" / "test_attestation.json"


def record_attestation(target_dir: Path | None = None, git_dir: Path | None = None) -> int:
    """Record a test pass receipt for current HEAD and tree."""
    state = get_current_state()
    if state is None:
        print("Warning: Could not determine Git state. Cannot record attestation.", file=sys.stderr)
        return 1

    head_sha, tree_sha, is_clean = state
    attestation_file = get_attestation_file(target_dir or git_dir)

    receipt = {
        "head_sha": head_sha,
        "tree_sha": tree_sha,
        "is_clean": is_clean,
        "timestamp": int(time.time()),
        "status": "PASSED",
    }

    try:
        attestation_file.parent.mkdir(parents=True, exist_ok=True)
        attestation_file.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
        print(f"✨ Recorded test attestation for commit {head_sha[:7]} (tree {tree_sha[:7]}).")
        return 0
    except OSError as exc:
        print(f"Warning: Failed to write attestation file: {exc}", file=sys.stderr)
        return 1


def verify_attestation(target_dir: Path | None = None, git_dir: Path | None = None, max_age: int = MAX_ATTESTATION_AGE_SECONDS) -> int:
    """Verify if a valid test pass receipt exists for current HEAD and tree. Return 0 on hit, 1 on miss."""
    attestation_file = get_attestation_file(target_dir or git_dir)
    if not attestation_file.exists():
        # Also check fallback locations if default wasn't found
        repo_root = Path(__file__).resolve().parent.parent
        candidates = [repo_root / ".private" / "test_attestation.json", repo_root / "build" / "test_attestation.json"]
        found = False
        for c in candidates:
            if c.exists():
                attestation_file = c
                found = True
                break
        if not found:
            return 1

    try:
        data = json.loads(attestation_file.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return 1

    state = get_current_state()
    if state is None:
        return 1

    head_sha, tree_sha, _ = state

    # Verify status
    if data.get("status") != "PASSED":
        return 1

    # Verify matching tree hash
    if data.get("tree_sha") != tree_sha:
        return 1

    # Verify matching commit SHA
    if data.get("head_sha") != head_sha:
        return 1

    # Verify freshness
    timestamp = data.get("timestamp", 0)
    age = int(time.time()) - timestamp
    if age > max_age or age < -60:  # Allow 60s clock skew
        return 1

    return 0


def invalidate_attestation(target_dir: Path | None = None, git_dir: Path | None = None) -> int:
    """Invalidate any existing test attestation receipt."""
    attestation_file = get_attestation_file(target_dir or git_dir)
    if attestation_file.exists():
        try:
            attestation_file.unlink()
        except OSError:
            pass
    # Also clean up secondary candidates
    repo_root = Path(__file__).resolve().parent.parent
    for d in [repo_root / ".private", repo_root / "build"]:
        candidate = d / "test_attestation.json"
        if candidate.exists():
            try:
                candidate.unlink()
            except OSError:
                pass
    return 0

def main() -> int:
    parser = argparse.ArgumentParser(description="Manage Git test attestation receipts.")
    parser.add_argument("action", choices=["record", "verify", "invalidate"], help="Action to perform")
    parser.add_argument("--max-age", type=int, default=MAX_ATTESTATION_AGE_SECONDS, help="Max valid receipt age in seconds")
    args = parser.parse_args()

    if args.action == "record":
        return record_attestation()
    elif args.action == "verify":
        return verify_attestation(max_age=args.max_age)
    elif args.action == "invalidate":
        return invalidate_attestation()

    return 1


if __name__ == "__main__":
    sys.exit(main())
