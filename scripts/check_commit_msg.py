#!/usr/bin/env python3
"""
check_commit_msg.py

Validates commit message subjects against Conventional Commits specification
and validates "Implements: REQ-..." trailers for requirements traceability.
Follows ASD-STE100 principles with deterministic diagnostics.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

CONVENTIONAL_PATTERN = re.compile(
    r"^(feat|fix|docs|style|refactor|perf|test|build|ci|chore|revert|release)"
    r"(\([a-z0-9_.,\-\/]+\))?!?: [^\s].+$"
)

REQ_ID_PATTERN = re.compile(r"^REQ-[A-Z0-9]{2,6}-\d{3,4}$")


def is_automated_or_special_commit(first_line: str) -> bool:
    """Return True if message matches standard Git merge/revert/fixup patterns."""
    prefixes = (
        "Merge ",
        "Revert ",
        "fixup! ",
        "squash! ",
        "amend! ",
        "Apply suggested changes",
    )
    return any(first_line.startswith(p) for p in prefixes)


def validate_commit_message(content: str) -> list[str]:
    """Validate commit message string and return list of diagnostic errors."""
    errors: list[str] = []
    lines = [line.strip() for line in content.splitlines() if not line.strip().startswith("#")]

    # Remove trailing/leading empty lines
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop()

    if not lines:
        errors.append("Commit message is empty.")
        return errors

    subject = lines[0]
    if is_automated_or_special_commit(subject):
        return []

    if not CONVENTIONAL_PATTERN.match(subject):
        errors.append(
            f"Subject '{subject}' does not follow Conventional Commits.\n"
            "  Expected format: <type>(<scope>)?: <description>\n"
            "  Allowed types:   feat, fix, docs, style, refactor, perf, test, build, ci, chore, revert, release\n"
            "  Example:         feat(core): implement BarrettGlaserAction (#14)"
        )

    # Validate trailers (e.g. Implements: REQ-ARCH-001)
    for line in lines[1:]:
        if line.startswith("Implements:"):
            raw_ids = line.split(":", 1)[1].strip()
            ids = [req_id.strip() for req_id in re.split(r"[,;\s]+", raw_ids) if req_id.strip()]
            for req_id in ids:
                if not REQ_ID_PATTERN.match(req_id):
                    errors.append(
                        f"Malformed requirement ID '{req_id}' in trailer '{line}'.\n"
                        "  Expected format: REQ-<SUBSYSTEM>-<NUMBER> (e.g. REQ-TRC-006)"
                    )

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Git commit message.")
    parser.add_argument("message_file", type=Path, help="Path to file containing commit message")
    args = parser.parse_args()

    if not args.message_file.exists():
        print(f"Error: Commit message file '{args.message_file}' not found.", file=sys.stderr)
        return 1

    content = args.message_file.read_text(encoding="utf-8")
    errors = validate_commit_message(content)

    if errors:
        print("\n❌ Invalid commit message:", file=sys.stderr)
        for err in errors:
            print(f"  • {err}", file=sys.stderr)
        print("", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
