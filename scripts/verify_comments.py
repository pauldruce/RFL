#!/usr/bin/env python3
"""
verify_comments.py - Enforce code comment standards and hygiene across RFL.

Rules enforced:
1. Unit test files (*tests*.cpp, test_*.py):
   - Must NOT contain Doxygen math/API tags (\\f$, \\f[, \\brief, @brief, @param, @return).
   - Multi-line block comments (/* ... */ or /** ... */) and line comments (//) are both permitted,
     provided they do not use Doxygen markup tags.
2. Header files (*.hpp):
   - @brief lines should not contain raw LaTeX formulas (\\f$ or \\f[).
     Keep brief summaries in clean text/code spans for IDE hover card readability;
     put complex LaTeX math in the body.
"""

import os
import re
import sys
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parent.parent

IGNORE_DIRS = {
    ".git",
    "build",
    "_deps",
    ".venv",
    "venv",
    "research",
    "src/legacy",
}

DOXYGEN_TAGS_IN_TESTS = [
    (re.compile(r"\\f[\$\[]"), "Doxygen LaTeX tag '\\f$' or '\\f[' (tests must use plain/Markdown math)"),
    (re.compile(r"[@\\]brief\b"), "Doxygen tag '@brief' / '\\brief' in unit test"),
    (re.compile(r"[@\\]param\b"), "Doxygen tag '@param' / '\\param' in unit test"),
    (re.compile(r"[@\\]return[s]?\b"), "Doxygen tag '@return' in unit test"),
]

BRIEF_WITH_LATEX = re.compile(r"[@\\]brief.*\\f[\$\[]")


def is_ignored(path: Path) -> bool:
    rel = path.relative_to(WORKSPACE_ROOT).as_posix()
    for ignored in IGNORE_DIRS:
        if rel == ignored or rel.startswith(f"{ignored}/"):
            return True
    return False


def verify_test_file(path: Path) -> list[str]:
    errors = []
    rel = path.relative_to(WORKSPACE_ROOT).as_posix()
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line_no, line in enumerate(f, start=1):
                for pattern, description in DOXYGEN_TAGS_IN_TESTS:
                    if pattern.search(line):
                        errors.append(
                            f"{rel}:{line_no} - Forbidden {description}\n    Line: {line.strip()}"
                        )
    except Exception as e:
        errors.append(f"{rel}: Failed to read file: {e}")
    return errors


def verify_header_file(path: Path) -> list[str]:
    errors = []
    rel = path.relative_to(WORKSPACE_ROOT).as_posix()
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line_no, line in enumerate(f, start=1):
                if BRIEF_WITH_LATEX.search(line):
                    errors.append(
                        f"{rel}:{line_no} - Raw LaTeX formula found in @brief line.\n"
                        f"    Keep @brief lines in plain text/code spans for IDE hover readability; put LaTeX math in the body.\n"
                        f"    Line: {line.strip()}"
                    )
    except Exception as e:
        errors.append(f"{rel}: Failed to read file: {e}")
    return errors


def main() -> int:
    all_errors = []
    files_checked = 0

    for root, dirs, files in os.walk(WORKSPACE_ROOT):
        dirs[:] = [d for d in dirs if not is_ignored(Path(root) / d)]

        for file in files:
            path = Path(root) / file
            if is_ignored(path):
                continue

            if path.suffix in {".cpp", ".cc"}:
                if "tests" in path.parts:
                    all_errors.extend(verify_test_file(path))
                    files_checked += 1
            elif path.suffix in {".hpp", ".h"}:
                if "tests" not in path.parts:
                    all_errors.extend(verify_header_file(path))
                    files_checked += 1

    if all_errors:
        print(f"❌ Comment verification failed with {len(all_errors)} issue(s):\n")
        for err in all_errors:
            print(f"  • {err}\n")
        print("See .agents/skills/technical-writing/SKILL.md#7-code-commenting-standards for guidelines.")
        return 1

    print(f"✅ Verified comment hygiene across {files_checked} C++ files. Zero violations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
