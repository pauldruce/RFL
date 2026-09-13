#!/usr/bin/env python3
"""Audits CHANGELOG.md to prevent releasing with missing PRs, placeholders, or stale tags."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CHANGELOG_PATH = REPO_ROOT / "CHANGELOG.md"


def get_latest_stable_tag() -> str | None:
    """Returns the most recent stable Git tag (vX.Y.Z)."""
    try:
        res = subprocess.run(
            ["git", "tag", "--sort=-v:refname", "--list", "v[0-9]*"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        tags = [t.strip() for t in res.stdout.splitlines() if t.strip()]
        stable_tags = [t for t in tags if re.match(r"^v[0-9]+\.[0-9]+\.[0-9]+$", t)]
        return stable_tags[0] if stable_tags else None
    except Exception:
        return None


def get_merged_prs_since(tag: str | None, excluded_prs: set[int] | None = None) -> dict[int, str]:
    """Finds all PR numbers and commit summaries merged since the given tag."""
    range_spec = f"{tag}..HEAD" if tag else "HEAD"
    prs: dict[int, str] = {}
    excluded = excluded_prs or set()
    try:
        res = subprocess.run(
            ["git", "log", range_spec, "--oneline"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        for line in res.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            # Look for (#123) pattern typical of squash merges; take the last match on the line
            matches = re.findall(r"\(#(\d+)\)", line)
            if matches:
                pr_num = int(matches[-1])
                if pr_num not in excluded:
                    prs[pr_num] = line
    except Exception as exc:
        print(f"⚠️ Warning: Could not extract git log: {exc}", file=sys.stderr)
    return prs


def audit_changelog(tag: str, previous_tag: str | None = None, excluded_prs: set[int] | None = None) -> list[str]:
    """Validates CHANGELOG.md against the specified tag and Git history."""
    errors: list[str] = []

    if not CHANGELOG_PATH.exists():
        return [f"CHANGELOG.md does not exist at {CHANGELOG_PATH}."]

    content = CHANGELOG_PATH.read_text(encoding="utf-8")

    # 1. Reject development versions
    if "devel" in tag.lower() or re.match(r"^v?0\.0\.0", tag):
        return [f"Tag '{tag}' is an unversioned development build. Official releases must use valid SemVer tags."]

    # 2. Extract base version (v0.2.0rc1 -> 0.2.0)
    semver_match = re.search(r"([0-9]+\.[0-9]+\.[0-9]+)", tag)
    if not semver_match:
        return [f"Tag '{tag}' does not contain a valid SemVer pattern (X.Y.Z)."]
    base_version = semver_match.group(1)

    # 3. Find target release section in CHANGELOG.md
    # Matches '## [v0.2.0]' or '## [0.2.0]'
    section_regex = re.compile(
        rf"^## \[(?:v)?{re.escape(base_version)}[^\]]*\].*?(?=^## |\Z)",
        re.MULTILINE | re.DOTALL,
    )
    section_match = section_regex.search(content)
    if not section_match:
        return [
            f"CHANGELOG.md has no entry for version '{base_version}'.\n"
            f"Run 'python3 scripts/draft_changelog.py --tag {tag} --write' and commit before tagging."
        ]

    section_text = section_match.group(0)

    # 4. Check for unresolved template placeholders
    placeholder_patterns = [
        re.compile(r"\bTODO\b", re.IGNORECASE),
        re.compile(r"\bTBD\b", re.IGNORECASE),
        re.compile(r"Add release highlights and breaking changes here", re.IGNORECASE),
        re.compile(r"\[Add.*?\]"),
    ]
    for pattern in placeholder_patterns:
        if pattern.search(section_text):
            errors.append(f"CHANGELOG.md section [{base_version}] contains unresolved placeholder: '{pattern.pattern}'")

    # 5. Check for missing PRs merged since previous tag
    effective_prev_tag = previous_tag or get_latest_stable_tag()
    # Don't compare against itself if the target tag is already the latest stable tag
    if effective_prev_tag == tag or effective_prev_tag == f"v{base_version}":
        effective_prev_tag = None

    merged_prs = get_merged_prs_since(effective_prev_tag, excluded_prs=excluded_prs)
    missing_prs: list[str] = []
    for pr_num, commit_line in merged_prs.items():
        # Check if #pr_num or /pull/pr_num appears in the section text
        pr_ref_pattern = rf"(#|/pull/){pr_num}\b"
        if not re.search(pr_ref_pattern, section_text):
            missing_prs.append(f"#{pr_num}: {commit_line}")

    if missing_prs:
        errors.append(
            f"CHANGELOG.md section [{base_version}] is missing {len(missing_prs)} merged PR(s) since {effective_prev_tag or 'initial commit'}:\n"
            + "\n".join(f"    - {item}" for item in missing_prs)
            + "\n    Tip: Self-reference all release and chore PRs in CHANGELOG.md (e.g. [#<pr>](https://github.com/pauldruce/RFL/pull/<pr>))."
        )

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit CHANGELOG.md for tag and PR completeness")
    parser.add_argument("--tag", type=str, required=True, help="Release tag to verify (e.g. v0.2.0 or v0.2.0rc1)")
    parser.add_argument("--previous-tag", type=str, default=None, help="Previous release tag (default: latest stable tag)")
    parser.add_argument(
        "--exclude-pr",
        type=int,
        action="append",
        default=[],
        help="PR number to exclude from verification (can be specified multiple times)",
    )

    args = parser.parse_args()

    errors = audit_changelog(args.tag, args.previous_tag, excluded_prs=set(args.exclude_pr))
    if errors:
        print(f"❌ CHANGELOG.md verification failed for tag '{args.tag}':", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print(f"✅ CHANGELOG.md verified successfully for release tag '{args.tag}'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
