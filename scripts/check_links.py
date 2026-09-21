#!/usr/bin/env python3
"""Validates internal relative links and GitHub Actions workflow badges in Markdown files."""

from __future__ import annotations

import re
import sys
from pathlib import Path

# Repository root directory
REPO_ROOT = Path(__file__).resolve().parent.parent

# Placeholders permitted in templates
PLACEHOLDER_PATTERN = re.compile(r"^(\.\.\.|\.\.\.|<.*>|TODO|TBD|FIXME|\.\.\./.*)$")


def check_markdown_file(file_path: Path) -> list[str]:
    """Checks a single Markdown file for broken relative links and stale workflow badges."""
    errors: list[str] = []
    try:
        content = file_path.read_text(encoding="utf-8")
    except Exception as exc:
        return [f"{file_path}: Failed to read file: {exc}"]

    # Match markdown links: [text](target) and image links: ![alt](target)
    link_pattern = re.compile(r"!?\[([^\]]*)\]\(([^)]+)\)")
    badge_pattern = re.compile(
        r"https://github\.com/[^/]+/[^/]+/actions/workflows/([^/]+)/badge\.svg"
    )
    workflow_run_pattern = re.compile(
        r"https://github\.com/[^/]+/[^/]+/actions/workflows/([^/)]+)"
    )

    is_template = "template" in file_path.name.lower()

    for line_no, line in enumerate(content.splitlines(), start=1):
        for match in link_pattern.finditer(line):
            target = match.group(2).strip()

            # Ignore templates containing placeholders
            if is_template and PLACEHOLDER_PATTERN.match(target):
                continue

            # 1. Check GitHub workflow badges
            badge_match = badge_pattern.search(target)
            if badge_match:
                workflow_name = badge_match.group(1)
                workflow_file = REPO_ROOT / ".github" / "workflows" / workflow_name
                if not workflow_file.exists():
                    errors.append(
                        f"{file_path}:{line_no}: Stale workflow badge references missing file: "
                        f"'.github/workflows/{workflow_name}'"
                    )
                continue

            # 2. Check workflow run links
            wf_match = workflow_run_pattern.search(target)
            if wf_match:
                workflow_name = wf_match.group(1).split("?")[0]
                workflow_file = REPO_ROOT / ".github" / "workflows" / workflow_name
                if not workflow_file.exists():
                    errors.append(
                        f"{file_path}:{line_no}: Stale workflow link references missing file: "
                        f"'.github/workflows/{workflow_name}'"
                    )
                continue

            # 3. Check relative local file links
            if (
                target.startswith("http://")
                or target.startswith("https://")
                or target.startswith("mailto:")
                or target.startswith("#")
            ):
                continue

            # Strip anchor (#section)
            file_target = target.split("#")[0]
            if not file_target:
                continue

            # Resolve path relative to containing directory
            resolved_target = (file_path.parent / file_target).resolve()

            # Ignore links pointing outside repository
            try:
                resolved_target.relative_to(REPO_ROOT)
            except ValueError:
                continue

            if not resolved_target.exists():
                errors.append(
                    f"{file_path}:{line_no}: Broken relative link to '{target}' "
                    f"(resolved to '{resolved_target}')"
                )

    return errors


def main() -> int:
    """Scans all Markdown files in the repository."""
    md_files = [
        REPO_ROOT / "README.md",
        REPO_ROOT / "CHANGELOG.md",
        REPO_ROOT / "ROADMAP.md",
        REPO_ROOT / "AGENTS.md",
    ]
    for path in (REPO_ROOT / "docs").rglob("*.md"):
        md_files.append(path)
    for path in (REPO_ROOT / "examples").rglob("*.md"):
        md_files.append(path)

    all_errors: list[str] = []
    for md_file in sorted(md_files):
        if md_file.exists():
            all_errors.extend(check_markdown_file(md_file))

    if all_errors:
        print(f"❌ Found {len(all_errors)} broken link(s) / stale badge(s):", file=sys.stderr)
        for err in all_errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print(f"✅ Verified {len(md_files)} markdown files. All internal links and workflow badges valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
