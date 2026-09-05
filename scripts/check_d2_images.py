#!/usr/bin/env python3
"""
check_d2_images.py

Linter and synchronisation tool for D2 diagram assets.
Verifies that all compiled SVG assets match their .d2 source specifications.
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_EXCLUDES = [
    "**/build/**",
    "**/research/**",
    "**/.worktrees/**",
    "**/venv/**",
    "**/.*/**",
]


def resolve_d2_binary(custom_path: str | None = None) -> Path:
    """Locate the d2 executable from argument, PATH, or mise shims."""
    if custom_path:
        p = Path(custom_path).expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise FileNotFoundError(f"Specified d2 binary not found or not executable: {custom_path}")

    # Check standard PATH
    path_d2 = shutil.which("d2")
    if path_d2:
        return Path(path_d2)

    # Check mise shim
    mise_d2 = Path.home() / ".local" / "share" / "mise" / "shims" / "d2"
    if mise_d2.is_file() and os.access(mise_d2, os.X_OK):
        return mise_d2

    # Check default install path
    default_d2 = Path("/usr/local/bin/d2")
    if default_d2.is_file() and os.access(default_d2, os.X_OK):
        return default_d2

    raise FileNotFoundError(
        "d2 executable not found on PATH or in ~/.local/share/mise/shims/d2.\n"
        "Install d2 using: mise use -g aqua:terrastruct/d2\n"
        "Or visit: https://d2lang.com/tour/install"
    )


def should_exclude(path: Path, excludes: list[str], root_dir: Path) -> bool:
    """Check if the given path matches any exclude glob pattern."""
    rel_str = str(path.relative_to(root_dir))
    for pat in excludes:
        if fnmatch.fnmatch(rel_str, pat) or fnmatch.fnmatch(f"**/{rel_str}", pat):
            return True
    return False


def find_d2_files(search_paths: list[str], excludes: list[str], root_dir: Path) -> list[Path]:
    """Find all relevant .d2 files according to search paths and exclusions."""
    d2_files: list[Path] = []

    if search_paths:
        for p_str in search_paths:
            p = Path(p_str).resolve()
            if p.is_file() and p.suffix == ".d2":
                d2_files.append(p)
            elif p.is_dir():
                for f in p.rglob("*.d2"):
                    d2_files.append(f)
    else:
        # Default: search docs/ and repository root
        for f in root_dir.rglob("*.d2"):
            if not should_exclude(f, excludes, root_dir):
                d2_files.append(f)

    return sorted(set(d2_files))


def compile_d2(d2_bin: Path, d2_file: Path) -> tuple[int, str, str]:
    """Compile a .d2 file to SVG format in-memory via stdout."""
    cmd = [str(d2_bin), str(d2_file), "-"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify or update compiled SVG diagrams generated from D2 sources."
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="Optional specific .d2 files or directories to check. Defaults to scanning the repository.",
    )
    parser.add_argument(
        "--fix",
        "-f",
        action="store_true",
        help="Automatically recompile and update missing or out-of-date SVG files.",
    )
    parser.add_argument(
        "--d2-path",
        default=None,
        help="Path to the d2 binary executable.",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help="Glob pattern to exclude (can be specified multiple times).",
    )

    args = parser.parse_args()
    root_dir = Path.cwd()

    try:
        d2_bin = resolve_d2_binary(args.d2_path)
    except FileNotFoundError as err:
        print(f"Error: {err}", file=sys.stderr)
        return 2

    excludes = DEFAULT_EXCLUDES + args.exclude
    d2_files = find_d2_files(args.files, excludes, root_dir)

    if not d2_files:
        print("No D2 files found to check.")
        return 0

    print(f"Checking {len(d2_files)} D2 diagram source(s) using {d2_bin}...")

    total = len(d2_files)
    up_to_date = 0
    stale = 0
    missing = 0
    errors = 0

    for d2_path in d2_files:
        rel_d2 = d2_path.relative_to(root_dir)
        svg_path = d2_path.with_suffix(".svg")
        rel_svg = svg_path.relative_to(root_dir)

        code, compiled_svg, err = compile_d2(d2_bin, d2_path)
        if code != 0:
            print(f"  [ERROR] {rel_d2}: Compilation failed:\n{err}", file=sys.stderr)
            errors += 1
            continue

        if not svg_path.exists():
            missing += 1
            if args.fix:
                svg_path.write_text(compiled_svg, encoding="utf-8")
                print(f"  [FIXED] {rel_svg} generated from {rel_d2}")
                up_to_date += 1
            else:
                print(f"  [MISSING] {rel_svg} does not exist for {rel_d2}")
            continue

        try:
            current_svg = svg_path.read_text(encoding="utf-8")
        except OSError as e:
            print(f"  [ERROR] Failed to read {rel_svg}: {e}", file=sys.stderr)
            errors += 1
            continue

        if current_svg == compiled_svg:
            print(f"  [OK] {rel_svg} is up to date with {rel_d2}")
            up_to_date += 1
        else:
            stale += 1
            if args.fix:
                svg_path.write_text(compiled_svg, encoding="utf-8")
                print(f"  [UPDATED] {rel_svg} recompiled from {rel_d2}")
                up_to_date += 1
            else:
                print(f"  [STALE] {rel_svg} is out of date with {rel_d2}")

    print("\n--- D2 Diagram Lint Summary ---")
    print(f"Total checked: {total}")
    print(f"Up to date:    {up_to_date}")
    print(f"Stale:         {stale}")
    print(f"Missing:       {missing}")
    print(f"Errors:        {errors}")

    if not args.fix and (stale > 0 or missing > 0 or errors > 0):
        print("\nFix required: Run 'python3 scripts/check_d2_images.py --fix' to regenerate stale/missing diagrams.")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
