#!/usr/bin/env python3
"""Drafts release notes and updates CHANGELOG.md locally using GitHub API or Git log."""

from __future__ import annotations

import argparse
import datetime
import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CHANGELOG_PATH = REPO_ROOT / "CHANGELOG.md"


def get_latest_git_tag(include_prereleases: bool = False) -> str | None:
    """Returns the most recent Git release tag (defaulting to stable tags)."""
    try:
        res = subprocess.run(
            ["git", "tag", "--sort=-v:refname", "--list", "v[0-9]*"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        tags = [t.strip() for t in res.stdout.splitlines() if t.strip()]
        if not include_prereleases:
            stable_tags = [t for t in tags if re.match(r"^v[0-9]+\.[0-9]+\.[0-9]+$", t)]
            if stable_tags:
                return stable_tags[0]
        return tags[0] if tags else None
    except Exception:
        return None


def fetch_github_generated_notes(tag: str, previous_tag: str | None) -> str | None:
    """Calls GitHub Releases generate-notes API via gh CLI."""
    cmd = [
        "gh",
        "api",
        "repos/pauldruce/RFL/releases/generate-notes",
        "-f",
        f"tag_name={tag}",
    ]
    if previous_tag:
        cmd.extend(["-f", f"previous_tag_name={previous_tag}"])

    try:
        res = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, check=True)
        data = json.loads(res.stdout)
        return data.get("body")
    except Exception:
        return None


def generate_notes_from_git_log(previous_tag: str | None) -> str:
    """Generates categorized release notes from git commits when offline."""
    range_spec = f"{previous_tag}..HEAD" if previous_tag else "HEAD"
    cmd = ["git", "log", range_spec, "--oneline", "--no-merges"]
    try:
        res = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, check=True)
        commits = [line.strip() for line in res.stdout.splitlines() if line.strip()]
    except Exception as exc:
        return f"* Failed to extract git log: {exc}"

    categories: dict[str, list[str]] = {
        "🚀 Features & Enhancements": [],
        "🐛 Bug Fixes": [],
        "🧰 Build & CI/CD Architecture": [],
        "📚 Documentation & Governance": [],
        "🔍 Other Changes": [],
    }

    for c in commits:
        parts = c.split(" ", 1)
        if len(parts) < 2:
            continue
        msg = parts[1]
        if msg.startswith("feat"):
            categories["🚀 Features & Enhancements"].append(f"* {msg}")
        elif msg.startswith("fix"):
            categories["🐛 Bug Fixes"].append(f"* {msg}")
        elif msg.startswith("ci") or msg.startswith("build"):
            categories["🧰 Build & CI/CD Architecture"].append(f"* {msg}")
        elif msg.startswith("docs"):
            categories["📚 Documentation & Governance"].append(f"* {msg}")
        elif not msg.startswith("chore"):
            categories["🔍 Other Changes"].append(f"* {msg}")

    lines: list[str] = []
    for title, items in categories.items():
        if items:
            lines.append(f"### {title}")
            lines.extend(items)
            lines.append("")

    return "\n".join(lines).strip()


def format_release_block(tag: str, notes: str, release_date: str) -> str:
    """Formats the release block according to Keep a Changelog."""
    # Ensure tag header has [vX.Y.Z] format
    tag_clean = tag if tag.startswith("v") else f"v{tag}"
    header = f"## [{tag_clean}] - {release_date}"

    return (
        f"{header}\n\n"
        f"### Highlights & Breaking Changes\n"
        f"* Add release highlights and breaking changes here.\n\n"
        f"{notes.strip()}\n"
    )


def update_changelog_file(release_block: str, tag: str) -> None:
    """Prepends or updates the release block in CHANGELOG.md."""
    tag_clean = tag if tag.startswith("v") else f"v{tag}"
    section_pattern = re.compile(rf"^## \[{re.escape(tag_clean)}\].*?(?=^## |\Z)", re.MULTILINE | re.DOTALL)

    initial_content = (
        "# Changelog\n\n"
        "All notable changes to the Random Fuzzy Library (RFL) are documented in this file.\n"
        "The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),\n"
        "and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).\n\n"
    )

    if not CHANGELOG_PATH.exists():
        new_content = initial_content + release_block + "\n"
        CHANGELOG_PATH.write_text(new_content, encoding="utf-8")
        print(f"✅ Created {CHANGELOG_PATH} with section [{tag_clean}].")
        return

    content = CHANGELOG_PATH.read_text(encoding="utf-8")

    # If this version already exists, replace it in-place
    if section_pattern.search(content):
        new_content = section_pattern.sub(release_block + "\n", content, count=1)
        CHANGELOG_PATH.write_text(new_content, encoding="utf-8")
        print(f"✅ Updated existing section [{tag_clean}] in {CHANGELOG_PATH}.")
        return

    # Otherwise prepend after the preamble
    match = re.search(r"^## ", content, re.MULTILINE)
    if match:
        idx = match.start()
        new_content = content[:idx] + release_block + "\n\n" + content[idx:]
    else:
        new_content = content.rstrip() + "\n\n" + release_block + "\n"

    CHANGELOG_PATH.write_text(new_content, encoding="utf-8")
    print(f"✅ Prepended section [{tag_clean}] to {CHANGELOG_PATH}.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Draft release notes and update CHANGELOG.md")
    parser.add_argument("--tag", type=str, default="v0.2.0", help="Target release tag (e.g. v0.2.0)")
    parser.add_argument("--previous-tag", type=str, default=None, help="Previous release tag (default: latest git tag)")
    parser.add_argument("--date", type=str, default=None, help="Release date (default: today)")
    parser.add_argument("--preview", action="store_true", help="Preview output without modifying CHANGELOG.md")
    parser.add_argument("--write", action="store_true", help="Write release notes to CHANGELOG.md")

    args = parser.parse_args()

    tag = args.tag
    previous_tag = args.previous_tag or get_latest_git_tag()
    release_date = args.date or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")

    # Fetch notes
    notes = fetch_github_generated_notes(tag, previous_tag)
    if not notes:
        print("ℹ️ GitHub API unavailable or returned empty. Falling back to local git log...", file=sys.stderr)
        notes = generate_notes_from_git_log(previous_tag)

    block = format_release_block(tag, notes, release_date)

    if args.preview or not args.write:
        print("=== RELEASE NOTES PREVIEW ===")
        print(block)
        print("=============================")
        if not args.write:
            print("\nRun with --write to write this block to CHANGELOG.md.")

    if args.write:
        update_changelog_file(block, tag)

    return 0


if __name__ == "__main__":
    sys.exit(main())
