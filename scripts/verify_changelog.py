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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tag_resolver import resolve_latest_stable_tag



def get_merged_prs_since(
    tag: str | None,
    end_ref: str = "HEAD",
    excluded_prs: set[int] | None = None,
    ignore_dependencies: bool = True,
) -> dict[int, str]:
    """Finds all PR numbers and commit summaries merged between tag and end_ref."""
    range_spec = f"{tag}..{end_ref}" if tag else end_ref
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
                if pr_num in excluded:
                    continue

                # Strip short commit hash from the beginning of oneline output
                summary = line.split(" ", 1)[1] if " " in line else line

                # Automated dependency updates do not require manual changelog entries
                if ignore_dependencies and summary.startswith(
                    ("ci(deps):", "build(deps):", "chore(deps):", "deps:")
                ):
                    continue

                prs[pr_num] = line
    except Exception as exc:
        print(f"⚠️ Warning: Could not extract git log: {exc}", file=sys.stderr)
    return prs


def audit_changelog(
    tag: str,
    previous_tag: str | None = None,
    until_ref: str | None = None,
    excluded_prs: set[int] | None = None,
    ignore_dependencies: bool = True,
) -> list[str]:
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
    effective_prev_tag = previous_tag or resolve_latest_stable_tag(excluding_tag=tag)

    if until_ref:
        effective_until = until_ref
    else:
        # If tag exists as a git ref (e.g. verifying an existing release), verify up to that tag.
        # Otherwise default to HEAD (e.g. verifying unreleased changes for a prospective tag).
        tag_exists = subprocess.run(
            ["git", "rev-parse", "--verify", "--quiet", tag],
            cwd=REPO_ROOT,
            capture_output=True,
        ).returncode == 0
        effective_until = tag if tag_exists else "HEAD"

    merged_prs = get_merged_prs_since(
        effective_prev_tag,
        end_ref=effective_until,
        excluded_prs=excluded_prs,
        ignore_dependencies=ignore_dependencies,
    )
    missing_prs: list[str] = []
    for pr_num, commit_line in merged_prs.items():
        # Check if #pr_num or /pull/pr_num appears in the section text
        pr_ref_pattern = rf"(#|/pull/){pr_num}\b"
        if not re.search(pr_ref_pattern, section_text):
            missing_prs.append(f"#{pr_num}: {commit_line}")

    if missing_prs:
        errors.append(
            f"CHANGELOG.md section [{base_version}] is missing {len(missing_prs)} merged PR(s) between {effective_prev_tag or 'initial commit'} and {effective_until}:\n"
            + "\n".join(f"    - {item}" for item in missing_prs)
            + "\n    Tip: Self-reference all release and chore PRs in CHANGELOG.md (e.g. [#<pr>](https://github.com/pauldruce/RFL/pull/<pr>))."
        )

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit CHANGELOG.md for tag and PR completeness")
    parser.add_argument("--tag", type=str, default=None, help="Release tag to verify (default: latest stable tag)")
    parser.add_argument("--previous-tag", type=str, default=None, help="Previous release tag (default: latest stable tag)")
    parser.add_argument(
        "--until",
        type=str,
        default=None,
        help="Git ref up to which commits should be verified (default: tag if it exists in git, otherwise HEAD)",
    )
    parser.add_argument(
        "--exclude-pr",
        type=int,
        action="append",
        default=[],
        help="PR number to exclude from verification (can be specified multiple times)",
    )
    parser.add_argument(
        "--include-deps",
        action="store_true",
        default=False,
        help="Require automated dependency updates (ci(deps), build(deps)) in CHANGELOG.md audit",
    )

    args = parser.parse_args()

    tag_to_verify = args.tag or resolve_latest_stable_tag()
    if not tag_to_verify:
        print("❌ Error: No Git tag specified and no stable release tags found.", file=sys.stderr)
        return 1

    errors = audit_changelog(
        tag_to_verify,
        args.previous_tag,
        until_ref=args.until,
        excluded_prs=set(args.exclude_pr),
        ignore_dependencies=not args.include_deps,
    )
    if errors:
        print(f"❌ CHANGELOG.md verification failed for tag '{tag_to_verify}':", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print(f"✅ CHANGELOG.md verified successfully for release tag '{tag_to_verify}'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
