#!/usr/bin/env python3
"""
test_test_attestation.py

Unit tests for test attestation receipt recording and verification.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

SCRIPT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import test_attestation


class TestTestAttestation(unittest.TestCase):
    """Test cryptographic tree-hash attestation mechanics."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.git_dir = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_verify_missing_file_returns_1(self) -> None:
        rc = test_attestation.verify_attestation(git_dir=self.git_dir)
        self.assertEqual(rc, 1)

    def test_verify_matching_state_returns_0(self) -> None:
        state = test_attestation.get_current_state()
        if state is None:
            self.skipTest("Not in a git repo")

        head_sha, tree_sha, is_clean = state
        attestation_file = self.git_dir / "test_attestation.json"
        attestation_file.write_text(
            json.dumps({
                "head_sha": head_sha,
                "tree_sha": tree_sha,
                "is_clean": is_clean,
                "timestamp": int(time.time()),
                "status": "PASSED",
            }),
            encoding="utf-8",
        )

        rc = test_attestation.verify_attestation(git_dir=self.git_dir)
        self.assertEqual(rc, 0)

    def test_verify_tree_mismatch_returns_1(self) -> None:
        state = test_attestation.get_current_state()
        if state is None:
            self.skipTest("Not in a git repo")

        head_sha, _, is_clean = state
        attestation_file = self.git_dir / "test_attestation.json"
        attestation_file.write_text(
            json.dumps({
                "head_sha": head_sha,
                "tree_sha": "different_tree_sha_12345",
                "is_clean": is_clean,
                "timestamp": int(time.time()),
                "status": "PASSED",
            }),
            encoding="utf-8",
        )

        rc = test_attestation.verify_attestation(git_dir=self.git_dir)
        self.assertEqual(rc, 1)

    def test_verify_expired_timestamp_returns_1(self) -> None:
        state = test_attestation.get_current_state()
        if state is None:
            self.skipTest("Not in a git repo")

        head_sha, tree_sha, is_clean = state
        attestation_file = self.git_dir / "test_attestation.json"
        attestation_file.write_text(
            json.dumps({
                "head_sha": head_sha,
                "tree_sha": tree_sha,
                "is_clean": is_clean,
                "timestamp": int(time.time()) - 10000,
                "status": "PASSED",
            }),
            encoding="utf-8",
        )

        rc = test_attestation.verify_attestation(git_dir=self.git_dir, max_age=3600)
        self.assertEqual(rc, 1)

    def test_verify_head_mismatch_returns_1(self) -> None:
        state = test_attestation.get_current_state()
        if state is None:
            self.skipTest("Not in a git repo")

        _, tree_sha, is_clean = state
        attestation_file = self.git_dir / "test_attestation.json"
        attestation_file.write_text(
            json.dumps({
                "head_sha": "different_head_sha_09876",
                "tree_sha": tree_sha,
                "is_clean": is_clean,
                "timestamp": int(time.time()),
                "status": "PASSED",
            }),
            encoding="utf-8",
        )

        rc = test_attestation.verify_attestation(git_dir=self.git_dir)
        self.assertEqual(rc, 1)

    def test_verify_status_not_passed_returns_1(self) -> None:
        state = test_attestation.get_current_state()
        if state is None:
            self.skipTest("Not in a git repo")

        head_sha, tree_sha, is_clean = state
        attestation_file = self.git_dir / "test_attestation.json"
        attestation_file.write_text(
            json.dumps({
                "head_sha": head_sha,
                "tree_sha": tree_sha,
                "is_clean": is_clean,
                "timestamp": int(time.time()),
                "status": "FAILED",
            }),
            encoding="utf-8",
        )

        rc = test_attestation.verify_attestation(git_dir=self.git_dir)
        self.assertEqual(rc, 1)

    def test_record_attestation_writes_valid_json(self) -> None:
        rc = test_attestation.record_attestation(target_dir=self.git_dir)
        self.assertEqual(rc, 0)
        attestation_file = self.git_dir / "test_attestation.json"
        self.assertTrue(attestation_file.exists())

        data = json.loads(attestation_file.read_text(encoding="utf-8"))
        self.assertEqual(data["status"], "PASSED")
        self.assertIn("head_sha", data)
        self.assertIn("tree_sha", data)

        # Now verify should pass
        verify_rc = test_attestation.verify_attestation(git_dir=self.git_dir)
        self.assertEqual(verify_rc, 0)

    def test_invalidate_removes_file(self) -> None:
        attestation_file = self.git_dir / "test_attestation.json"
        attestation_file.write_text("{}", encoding="utf-8")
        self.assertTrue(attestation_file.exists())

        test_attestation.invalidate_attestation(git_dir=self.git_dir)
        self.assertFalse(attestation_file.exists())


if __name__ == "__main__":
    unittest.main()
