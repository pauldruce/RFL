#!/usr/bin/env python3
"""
tag_resolver.py

Single-purpose module for resolving and inspecting Git release tags and
semantic versions across RFL automation scripts.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Repository root directory
REPO_ROOT = Path(__file__).resolve().parent.parent

# Regex matching standard stable semantic version release tags (e.g. v0.2.0, v1.0.4)
STABLE_TAG_PATTERN = re.compile(r"^v?([0-9]+)\.([0-9]+)\.([0-9]+)$")

# Regex matching semantic version release tags including pre-releases (e.g. v0.2.0rc1, v1.0.0-beta.1)
SEMVER_TAG_PATTERN = re.compile(r"^v?([0-9]+)\.([0-9]+)\.([0-9]+)(?:[-.]?([a-zA-Z0-9_.-]+))?$")


def is_stable_tag(tag: str) -> bool:
    """Return True if the tag matches a stable release format (vMAJOR.MINOR.PATCH)."""
    return bool(STABLE_TAG_PATTERN.match(tag.strip()))


def is_prerelease_tag(tag: str) -> bool:
    """Return True if the tag matches a valid semver tag containing pre-release identifiers."""
    match = SEMVER_TAG_PATTERN.match(tag.strip())
    if not match:
        return False
    prerelease_part = match.group(4)
    return bool(prerelease_part)


def parse_semver(tag: str) -> tuple[int, int, int, str | None]:
    """Parse a Git tag into (major, minor, patch, prerelease_suffix).

    Raises:
        ValueError: If the tag does not conform to semantic versioning.
    """
    clean = tag.strip()
    match = SEMVER_TAG_PATTERN.match(clean)
    if not match:
        raise ValueError(f"Tag '{tag}' is not a valid semantic version.")
    major = int(match.group(1))
    minor = int(match.group(2))
    patch = int(match.group(3))
    prerelease = match.group(4)
    return major, minor, patch, prerelease


def list_release_tags(
    repo_root: Path | None = None,
    include_prereleases: bool = True,
) -> list[str]:
    """Return sorted list of release tags from Git history, descending by version."""
    root = repo_root or REPO_ROOT
    try:
        res = subprocess.run(
            ["git", "tag", "--sort=-v:refname", "--list", "v[0-9]*"],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        )
        raw_tags = [t.strip() for t in res.stdout.splitlines() if t.strip()]
        if include_prereleases:
            return [t for t in raw_tags if SEMVER_TAG_PATTERN.match(t)]
        return [t for t in raw_tags if is_stable_tag(t)]
    except Exception:
        return []


def resolve_latest_tag(
    repo_root: Path | None = None,
    include_prereleases: bool = True,
    excluding_tag: str | None = None,
) -> str | None:
    """Return the most recent Git release tag (including pre-releases if requested)."""
    tags = list_release_tags(repo_root=repo_root, include_prereleases=include_prereleases)
    if excluding_tag:
        clean_ex = excluding_tag.strip()
        tags = [t for t in tags if t != clean_ex]
    return tags[0] if tags else None


def resolve_latest_stable_tag(
    repo_root: Path | None = None,
    excluding_tag: str | None = None,
) -> str | None:
    """Return the most recent stable Git release tag (vX.Y.Z), optionally excluding a tag."""
    tags = list_release_tags(repo_root=repo_root, include_prereleases=False)
    if excluding_tag:
        clean_ex = excluding_tag.strip()
        base_match = STABLE_TAG_PATTERN.match(clean_ex.split("rc")[0].split("-")[0])
        base_version = f"v{base_match.group(1)}.{base_match.group(2)}.{base_match.group(3)}" if base_match else None
        tags = [t for t in tags if t != clean_ex and t != base_version]
    return tags[0] if tags else None


def main() -> int:
    """CLI entry point for querying release tags."""
    parser = argparse.ArgumentParser(description="Query and resolve RFL Git release tags")
    parser.add_argument(
        "--stable",
        action="store_true",
        help="Resolve latest stable tag (excluding pre-releases)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="List all release tags in descending version order",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        default=None,
        help="Tag to exclude from resolution",
    )

    args = parser.parse_args()

    if args.all:
        tags = list_release_tags(include_prereleases=not args.stable)
        for t in tags:
            print(t)
        return 0

    if args.stable:
        tag = resolve_latest_stable_tag(excluding_tag=args.exclude)
    else:
        tag = resolve_latest_tag(include_prereleases=True, excluding_tag=args.exclude)

    if tag:
        print(tag)
        return 0

    print("No release tags found.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
