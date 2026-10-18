#!/usr/bin/env python3
"""
Unit tests for scripts/verify_comments.py (Code Comment Standards and Hygiene).

Verifies:
- Detection of forbidden Doxygen tags in unit test files (\\f$, \\f[, @brief, @param, @return)
- Detection of raw LaTeX formulas in header @brief lines
- Preservation of clean Markdown/plain-text comments in unit tests and headers
- Directory exclusion logic for ignored paths (build, research, legacy, etc.)
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import verify_comments  # noqa: E402


class TestVerifyComments(unittest.TestCase):
    """Unit tests for code comment standards and hygiene enforcement."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_root = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_clean_unit_test_passes(self) -> None:
        test_file = self.test_root / "tDelta.cpp"
        test_file.write_text(
            "// Normal comment explaining why tolerance is 1e-12\n"
            "/* Block comment with ASCII math: Tr(D^2) */\n"
            "TEST(DeltaTests, MatchesTrace) {\n"
            "    EXPECT_TRUE(true);\n"
            "}\n",
            encoding="utf-8",
        )
        errors = verify_comments.verify_test_file(test_file)
        self.assertEqual(errors, [])

    def test_forbidden_doxygen_in_unit_test_detected(self) -> None:
        forbidden_cases = [
            ("// \\f$ \\Delta S \\f$", "LaTeX tag"),
            ("// \\f[ \\Delta S \\f]", "LaTeX tag"),
            ("/// @brief Test description", "@brief"),
            ("/// \\brief Test description", "\\brief"),
            ("/// @param x Parameter", "@param"),
            ("/// @return Result", "@return"),
        ]
        for comment_line, label in forbidden_cases:
            test_file = self.test_root / "test_case.cpp"
            test_file.write_text(
                f"{comment_line}\n"
                "TEST(Suite, Case) {}\n",
                encoding="utf-8",
            )
            errors = verify_comments.verify_test_file(test_file)
            self.assertTrue(
                len(errors) > 0,
                f"Expected violation for '{comment_line}' ({label}), but got none.",
            )

    def test_clean_header_passes(self) -> None:
        header_file = self.test_root / "Operator.hpp"
        header_file.write_text(
            "/**\n"
            " * @brief Calculates the spectral action gradient in ASCII/code spans.\n"
            " *\n"
            " * The exact formula is given by:\n"
            " * \\f[ \\nabla S = 2(D + D^\\dagger) \\f]\n"
            " */\n"
            "class Operator {};\n",
            encoding="utf-8",
        )
        errors = verify_comments.verify_header_file(header_file)
        self.assertEqual(errors, [])

    def test_latex_in_header_brief_detected(self) -> None:
        header_file = self.test_root / "BadHeader.hpp"
        header_file.write_text(
            "/**\n"
            " * @brief Calculates \\f$ \\mathrm{Tr}(D^2) \\f$ directly.\n"
            " */\n"
            "class BadHeader {};\n",
            encoding="utf-8",
        )
        errors = verify_comments.verify_header_file(header_file)
        self.assertTrue(
            len(errors) > 0,
            "Expected violation for LaTeX math in @brief line, but got none.",
        )

    def test_ignored_directory_filtering(self) -> None:
        repo_root = verify_comments.REPO_ROOT
        self.assertTrue(verify_comments.is_ignored(repo_root / "build"))
        self.assertTrue(verify_comments.is_ignored(repo_root / "research"))
        self.assertTrue(verify_comments.is_ignored(repo_root / ".git"))
        self.assertTrue(verify_comments.is_ignored(repo_root / "src/legacy"))
        self.assertFalse(verify_comments.is_ignored(repo_root / "src/core"))


if __name__ == "__main__":
    unittest.main()
