#!/usr/bin/env python3
"""Validates code comment standards and hygiene across C++ source and header files.

Enforces rules defined in .agents/skills/technical-writing/SKILL.md:
1. Unit test files (*tests* or test files):
   - Must NOT contain Doxygen math/API tags (\\f$, \\f[, \\brief, @brief, @param, @return).
   - Standard line (//) and block (/* ... */) comments with ASCII/Markdown math are permitted.
2. Public header files (*.hpp, *.h):
   - @brief lines must not contain raw LaTeX formulas (\\f$ or \\f[).
   - Keep brief summaries in clean text/code spans for IDE hover readability; put LaTeX math in the body.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

# Repository root directory
REPO_ROOT = Path(__file__).resolve().parent.parent

# Directories ignored by the scanner
IGNORE_DIRS = {
    ".git",
    ".worktrees",
    "build",
    "_deps",
    ".venv",
    "venv",
    "research",
    "src/legacy",
}

# Forbidden Doxygen tags in unit test files
DOXYGEN_TAGS_IN_TESTS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\\f[\$\[]"), "Doxygen LaTeX tag '\\f$' or '\\f[' (tests must use plain/Markdown math)"),
    (re.compile(r"[@\\]brief\b"), "Doxygen tag '@brief' / '\\brief' in unit test"),
    (re.compile(r"[@\\]param\b"), "Doxygen tag '@param' / '\\param' in unit test"),
    (re.compile(r"[@\\]return[s]?\b"), "Doxygen tag '@return' in unit test"),
]

# Detect LaTeX formulas in @brief lines of header files
BRIEF_WITH_LATEX = re.compile(r"[@\\]brief.*\\f[\$\[]")


def is_ignored(path: Path) -> bool:
    """Checks if a given path matches any ignored directory."""
    try:
        rel = path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return True
    for ignored in IGNORE_DIRS:
        if rel == ignored or rel.startswith(f"{ignored}/"):
            return True
    return False


def verify_test_file(path: Path) -> list[str]:
    """Validates comment standards for a unit test file."""
    errors: list[str] = []
    try:
        rel = path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        rel = path.as_posix()
    try:
        content = path.read_text(encoding="utf-8")
    except Exception as exc:
        return [f"{rel}: Failed to read file: {exc}"]

    for line_no, line in enumerate(content.splitlines(), start=1):
        for pattern, description in DOXYGEN_TAGS_IN_TESTS:
            if pattern.search(line):
                errors.append(
                    f"{rel}:{line_no}: Forbidden {description}: '{line.strip()}'"
                )
    return errors


def verify_header_file(path: Path) -> list[str]:
    """Validates comment standards for a public header file."""
    errors: list[str] = []
    try:
        rel = path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        rel = path.as_posix()
    try:
        content = path.read_text(encoding="utf-8")
    except Exception as exc:
        return [f"{rel}: Failed to read file: {exc}"]

    for line_no, line in enumerate(content.splitlines(), start=1):
        if BRIEF_WITH_LATEX.search(line):
            errors.append(
                f"{rel}:{line_no}: Raw LaTeX formula in @brief line (keep @brief in plain text for IDE hover): '{line.strip()}'"
            )
    return errors


def main() -> int:
    """Scans C++ files and validates comment hygiene."""
    parser = argparse.ArgumentParser(
        description="Verify code comment standards and hygiene across C++ files."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        default=[REPO_ROOT / "src"],
        help="Files or directories to scan (defaults to src/)",
    )
    args = parser.parse_args()

    files_to_check: list[Path] = []
    for path in args.paths:
        target = path if path.is_absolute() else (REPO_ROOT / path)
        if target.is_file() and target.suffix in {".cpp", ".hpp", ".h", ".cc"}:
            if not is_ignored(target):
                files_to_check.append(target)
        elif target.is_dir():
            for root, dirs, files in os.walk(target):
                dirs[:] = [d for d in dirs if not is_ignored(Path(root) / d)]
                for file in files:
                    file_path = Path(root) / file
                    if file_path.suffix in {".cpp", ".hpp", ".h", ".cc"} and not is_ignored(file_path):
                        files_to_check.append(file_path)

    files_to_check.sort()
    all_errors: list[str] = []
    files_checked = 0

    for file_path in files_to_check:
        if file_path.suffix in {".cpp", ".cc"}:
            if "tests" in file_path.parts:
                all_errors.extend(verify_test_file(file_path))
                files_checked += 1
        elif file_path.suffix in {".hpp", ".h"}:
            if "tests" not in file_path.parts:
                all_errors.extend(verify_header_file(file_path))
                files_checked += 1

    if all_errors:
        print(f"❌ Found {len(all_errors)} comment violation(s):", file=sys.stderr)
        for err in all_errors:
            print(f"  • {err}", file=sys.stderr)
        return 1

    print(f"✅ Verified comment hygiene across {files_checked} C++ files. Zero violations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
