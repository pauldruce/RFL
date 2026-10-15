#!/usr/bin/env python3
"""scripts/tests/test_verify_tdd_commits.py

Unit test suite for the atomic Red-Green TDD commit protocol inspector.
Verifies Red Anchor resolution (parent of first feat: commit), handling of
multiple test: commits, chronological inversion prevention, and dry-run mode.

Conforms to EP-5 (REQ-TRC-007).
Uses standard Python library only (runnable via unittest or pytest).
"""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

# Add scripts directory to module search path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
import sys
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import verify_tdd_commits as vtc


class TestRedGreenAnchorResolution(unittest.TestCase):
    """Test deterministic Red Anchor resolution under EP-5 (REQ-TRC-007)."""

    def _init_git_repo(self, tmp_path: Path) -> None:
        """Initialise a temporary git repository with initial commit."""
        subprocess.run(["git", "init", "-b", "main"], cwd=tmp_path, check=True, capture_output=True)
        subprocess.run(["git", "config", "user.name", "Tester"], cwd=tmp_path, check=True)
        subprocess.run(["git", "config", "user.email", "tester@example.com"], cwd=tmp_path, check=True)

        readme = tmp_path / "README.md"
        readme.write_text("# Test Repo\n", encoding="utf-8")
        subprocess.run(["git", "add", "README.md"], cwd=tmp_path, check=True)
        subprocess.run(["git", "commit", "-m", "chore: initial commit"], cwd=tmp_path, check=True, capture_output=True)

    def test_two_commit_canonical_sequence(self) -> None:
        """Verify standard two-commit PR: Commit 1 (test:), Commit 2 (feat:)."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            self._init_git_repo(tmp_root)

            # Create feature branch
            subprocess.run(["git", "checkout", "-b", "feature"], cwd=tmp_root, check=True, capture_output=True)

            # Commit 1: test: add stubs and unit test
            test_file = tmp_root / "test_kernel.cpp"
            test_file.write_text("// stub test\n", encoding="utf-8")
            subprocess.run(["git", "add", "test_kernel.cpp"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "test: add stubs and unit test"], cwd=tmp_root, check=True, capture_output=True)
            test_commit_sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_root, capture_output=True, text=True).stdout.strip()

            # Commit 2: feat: implement math kernel
            impl_file = tmp_root / "kernel.cpp"
            impl_file.write_text("// kernel implementation\n", encoding="utf-8")
            subprocess.run(["git", "add", "kernel.cpp"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "feat: implement math kernel\n\nImplements: REQ-PHY-001"], cwd=tmp_root, check=True, capture_output=True)
            feat_commit_sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_root, capture_output=True, text=True).stdout.strip()

            red_sha, green_sha, errors = vtc.resolve_red_and_green_anchors("main", tmp_root)
            self.assertEqual(len(errors), 0, f"Errors: {errors}")
            self.assertEqual(red_sha, test_commit_sha, "Red Anchor must be immediate parent of first feat:")
            self.assertEqual(green_sha, feat_commit_sha, "Green Anchor must be HEAD")

    def test_chronological_inversion_prevention(self) -> None:
        """Verify parent-of-first-feat resolves correctly when test: commit follows feat: commit."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            self._init_git_repo(tmp_root)

            subprocess.run(["git", "checkout", "-b", "feature"], cwd=tmp_root, check=True, capture_output=True)

            # Commit 1: test: initial stub
            f1 = tmp_root / "f1.txt"
            f1.write_text("1\n", encoding="utf-8")
            subprocess.run(["git", "add", "f1.txt"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "test: initial stubs"], cwd=tmp_root, check=True, capture_output=True)

            # Commit 2: test: fix syntax error
            f2 = tmp_root / "f2.txt"
            f2.write_text("2\n", encoding="utf-8")
            subprocess.run(["git", "add", "f2.txt"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "test: fix syntax error and adjust stubs"], cwd=tmp_root, check=True, capture_output=True)
            true_red_sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_root, capture_output=True, text=True).stdout.strip()

            # Commit 3: feat: implement kernel
            f3 = tmp_root / "f3.txt"
            f3.write_text("3\n", encoding="utf-8")
            subprocess.run(["git", "add", "f3.txt"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "feat: implement mathematical kernel"], cwd=tmp_root, check=True, capture_output=True)

            # Commit 4: test: fix warning (occurs AFTER feat: commit)
            f4 = tmp_root / "f4.txt"
            f4.write_text("4\n", encoding="utf-8")
            subprocess.run(["git", "add", "f4.txt"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "test: fix compiler warning in test suite"], cwd=tmp_root, check=True, capture_output=True)
            head_sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_root, capture_output=True, text=True).stdout.strip()

            red_sha, green_sha, errors = vtc.resolve_red_and_green_anchors("main", tmp_root)
            self.assertEqual(len(errors), 0, f"Errors: {errors}")
            # Resolution 2 selects true_red_sha (parent of first feat:), NOT commit 4!
            self.assertEqual(red_sha, true_red_sha, "Must select parent of first feat:, avoiding chronological inversion")
            self.assertEqual(green_sha, head_sha)

    def test_refactoring_pr_without_feat_commit(self) -> None:
        """Verify that PR without feat: commit returns None (treated as non-blocking refactoring)."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            self._init_git_repo(tmp_root)

            subprocess.run(["git", "checkout", "-b", "chore-branch"], cwd=tmp_root, check=True, capture_output=True)

            f1 = tmp_root / "cleanup.txt"
            f1.write_text("clean\n", encoding="utf-8")
            subprocess.run(["git", "add", "cleanup.txt"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "refactor: clean up comments"], cwd=tmp_root, check=True, capture_output=True)

            red_sha, green_sha, errors = vtc.resolve_red_and_green_anchors("main", tmp_root)
            self.assertEqual(len(errors), 0)
            self.assertIsNone(red_sha)
            self.assertIsNone(green_sha)

            # Running verify_tdd_protocol in dry run returns 0 (non-blocking pass)
            ret = vtc.verify_tdd_protocol("main", tmp_root, dry_run=True)
            self.assertEqual(ret, 0)

    def test_dry_run_mode(self) -> None:
        """Verify dry-run mode returns 0 when valid Red/Green commits exist."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            self._init_git_repo(tmp_root)

            subprocess.run(["git", "checkout", "-b", "feature"], cwd=tmp_root, check=True, capture_output=True)

            f1 = tmp_root / "t.cpp"
            f1.write_text("test\n", encoding="utf-8")
            subprocess.run(["git", "add", "t.cpp"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "test: add stub test"], cwd=tmp_root, check=True, capture_output=True)

            f2 = tmp_root / "k.cpp"
            f2.write_text("kernel\n", encoding="utf-8")
            subprocess.run(["git", "add", "k.cpp"], cwd=tmp_root, check=True)
            subprocess.run(["git", "commit", "-m", "feat: add implementation"], cwd=tmp_root, check=True, capture_output=True)

            ret = vtc.verify_tdd_protocol("main", tmp_root, dry_run=True)
            self.assertEqual(ret, 0)


if __name__ == "__main__":
    unittest.main()
