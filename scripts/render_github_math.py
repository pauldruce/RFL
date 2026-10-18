#!/usr/bin/env python3
"""Validates Markdown mathematical rendering against the live GitHub Markdown API.

Requires network access (BypassSandbox: true if executed within an agent sandbox).
Rate limit: 60 requests/hour unauthenticated. Set GITHUB_TOKEN environment variable for higher limits.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def check_github_rendered_math(file_path: Path) -> list[str]:
    """Renders a Markdown file via GitHub's API and detects broken/unrendered LaTeX math."""
    try:
        text = file_path.read_text(encoding="utf-8")
    except Exception as exc:
        return [f"{file_path}: Failed to read file: {exc}"]

    url = "https://api.github.com/markdown"
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "RFL-Math-Validator",
    }
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"

    payload = json.dumps({
        "text": text,
        "mode": "gfm",
        "context": "pauldruce/RFL",
    }).encode("utf-8")

    req = urllib.request.Request(url, data=payload, headers=headers)

    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            html = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        return [f"{file_path}: GitHub API error HTTP {exc.code}: {exc.reason}"]
    except Exception as exc:
        return [f"{file_path}: Network error while calling GitHub API: {exc}"]

    errors: list[str] = []

    # Find all paragraph, heading, list item, and table cell blocks
    blocks = re.findall(r"<(?:p|h[1-6]|li|td|th)[^>]*>.*?</(?:p|h[1-6]|li|td|th)>", html, re.DOTALL)

    for block in blocks:
        # Strip rendered math tags (<math-renderer ...>...</math-renderer>)
        stripped = re.sub(r"<math-renderer.*?</math-renderer>", "[MATH]", block, flags=re.DOTALL)
        # Strip code elements (<code>...</code>, <pre>...</pre>)
        stripped = re.sub(r"<code>.*?</code>", "[CODE]", stripped, flags=re.DOTALL)

        # Check for unrendered dollar signs (excluding escaped \$)
        dollars = re.findall(r"(?<!\\)\$", stripped)
        if dollars:
            clean_snippet = re.sub(r"\s+", " ", stripped).strip()[:140]
            errors.append(f"{file_path}: Unrendered '$' in rendered HTML block: '{clean_snippet}'")
            continue

        # Check for LaTeX syntax accidentally eaten by Markdown italics (<em>...</em>)
        em_matches = re.findall(r"<em>(.*?)</em>", stripped)
        for em_content in em_matches:
            if any(char in em_content for char in ["\\", "{", "}", "^", "_"]):
                clean_snippet = re.sub(r"\s+", " ", stripped).strip()[:140]
                errors.append(
                    f"{file_path}: LaTeX subscript mangled into Markdown italics '<em>{em_content}</em>': '{clean_snippet}'"
                )
                break

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify math rendering in Markdown files via GitHub's API."
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="Markdown files to validate via GitHub API",
    )
    args = parser.parse_args()

    all_errors: list[str] = []
    checked_count = 0

    for path in args.paths:
        if path.is_file() and path.suffix == ".md":
            checked_count += 1
            print(f"Validating {path} via GitHub Markdown API...", file=sys.stderr)
            errs = check_github_rendered_math(path)
            all_errors.extend(errs)

    if all_errors:
        print(f"\n❌ Found {len(all_errors)} math rendering violation(s) on GitHub:", file=sys.stderr)
        for err in all_errors:
            print(f"  • {err}", file=sys.stderr)
        return 1

    print(f"\n✅ All math rendered cleanly on GitHub across {checked_count} file(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
