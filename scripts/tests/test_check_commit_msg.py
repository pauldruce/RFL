#!/usr/bin/env python3
"""
test_check_commit_msg.py

Unit tests for Conventional Commits and REQ trailer validation.
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

SCRIPT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import check_commit_msg


class TestCheckCommitMsg(unittest.TestCase):
    """Test Conventional Commits and requirement trailer parsing."""

    def test_valid_subjects(self) -> None:
        valid_messages = [
            "feat(core): implement DiracOperator value semantics",
            "fix: correct sign convention in Clifford algebra traces",
            "docs(changelog): record PR #84 under v0.3.0 features",
            "test(physics): add Barrett-Glaser action unit test",
            "refactor(build): adopt canonical include directory layout",
            "chore: update dependencies",
            "ci(plan): optimise target impact analysis",
            "feat!: break legacy API for non-blocking stepper",
            "feat(action)!: rewrite action interface",
            "Merge branch 'main' into feat/issue-14",
            "Revert \"feat(core): temporary change\"",
        ]
        for msg in valid_messages:
            errors = check_commit_msg.validate_commit_message(msg)
            self.assertEqual(errors, [], f"Expected '{msg}' to be valid")

    def test_invalid_subjects(self) -> None:
        invalid_messages = [
            "fixed stuff",
            "WIP: testing",
            "feat:",
            "feat(core):",
            "Random commit message without type prefix",
            "FEAT(core): uppercase is not allowed",
            "docs: ",
        ]
        for msg in invalid_messages:
            errors = check_commit_msg.validate_commit_message(msg)
            self.assertTrue(len(errors) > 0, f"Expected '{msg}' to be invalid")

    def test_valid_trailers(self) -> None:
        msg = (
            "feat(traceability): add shift-left git hook validation\n\n"
            "Detailed explanation of hook infrastructure.\n\n"
            "Implements: REQ-TRC-006\n"
        )
        errors = check_commit_msg.validate_commit_message(msg)
        self.assertEqual(errors, [])

    def test_multiple_valid_trailers(self) -> None:
        msg = (
            "feat(traceability): add multiple requirements\n\n"
            "Implements: REQ-TRC-006, REQ-CI-001\n"
        )
        errors = check_commit_msg.validate_commit_message(msg)
        self.assertEqual(errors, [])

    def test_invalid_trailers(self) -> None:
        msg = (
            "feat(traceability): add requirement\n\n"
            "Implements: INVALID_REQ\n"
        )
        errors = check_commit_msg.validate_commit_message(msg)
        self.assertTrue(any("Malformed requirement ID" in e for e in errors))

    def test_empty_message(self) -> None:
        errors = check_commit_msg.validate_commit_message("")
        self.assertTrue(any("empty" in e for e in errors))


if __name__ == "__main__":
    unittest.main()
