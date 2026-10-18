#!/usr/bin/env python3
"""
check_configs.py

Validates repository configuration files (TOML, JSON, YAML) for syntax integrity.
Follows ASD-STE100 principles with deterministic compiler-grade error reporting.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

try:
    import tomllib
except ImportError:
    import tomli as tomllib  # type: ignore

try:
    import yaml
except ImportError:
    import glob
    repo_root = Path(__file__).resolve().parent.parent
    for p in glob.glob(str(repo_root / ".venv/lib/python*/site-packages")):
        if p not in sys.path:
            sys.path.insert(0, p)
    try:
        import yaml
    except ImportError:
        yaml = None


IGNORE_DIRS = {
    ".git",
    ".venv",
    "venv",
    "build",
    "cmake-build",
    "dist",
    "_skbuild",
    ".skbuild",
    ".worktrees",
    ".private",
    "node_modules",
}


def should_skip(path: Path) -> bool:
    """Return True if path contains any ignored directory components."""
    return any(part in IGNORE_DIRS for part in path.parts)


def check_toml_file(path: Path) -> list[str]:
    """Validate a single TOML file."""
    errors: list[str] = []
    try:
        with open(path, "rb") as f:
            tomllib.load(f)
    except Exception as exc:
        errors.append(f"{path}:1: [TOML] Syntax error: {exc}")
    return errors


def check_json_file(path: Path) -> list[str]:
    """Validate a single JSON file."""
    errors: list[str] = []
    try:
        with open(path, "r", encoding="utf-8") as f:
            json.load(f)
    except Exception as exc:
        errors.append(f"{path}:1: [JSON] Syntax error: {exc}")
    return errors


def check_yaml_file(path: Path) -> list[str]:
    """Validate a single YAML file."""
    errors: list[str] = []
    if yaml is None:
        errors.append(f"{path}:1: [YAML] PyYAML is not installed; cannot validate YAML syntax.")
        return errors
    try:
        with open(path, "r", encoding="utf-8") as f:
            yaml.safe_load(f)
    except Exception as exc:
        errors.append(f"{path}:1: [YAML] Syntax error: {exc}")
    return errors


def find_files(repo_root: Path) -> tuple[list[Path], list[Path], list[Path]]:
    """Discover all non-ignored TOML, JSON, and YAML files."""
    toml_files: list[Path] = []
    json_files: list[Path] = []
    yaml_files: list[Path] = []

    for root, dirs, files in os.walk(repo_root):
        dirs[:] = [d for d in dirs if d not in IGNORE_DIRS]
        for f in files:
            p = Path(root) / f
            if should_skip(p):
                continue
            suffix = p.suffix.lower()
            if suffix == ".toml":
                toml_files.append(p)
            elif suffix == ".json":
                json_files.append(p)
            elif suffix in {".yml", ".yaml"}:
                yaml_files.append(p)

    return sorted(toml_files), sorted(json_files), sorted(yaml_files)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate TOML, JSON, and YAML config syntax.")
    parser.add_argument("--root", type=Path, default=Path("."), help="Repository root directory")
    args = parser.parse_args()

    repo_root = args.root.resolve()
    toml_files, json_files, yaml_files = find_files(repo_root)

    total_files = len(toml_files) + len(json_files) + len(yaml_files)
    all_errors: list[str] = []

    for f in toml_files:
        all_errors.extend(check_toml_file(f))

    for f in json_files:
        all_errors.extend(check_json_file(f))

    for f in yaml_files:
        all_errors.extend(check_yaml_file(f))

    if all_errors:
        for err in all_errors:
            print(err, file=sys.stderr)
        print(f"\n❌ Configuration audit failed: {len(all_errors)} errors in {total_files} files.", file=sys.stderr)
        return 1

    print(f"✅ Verified {total_files} configuration files ({len(toml_files)} TOML, {len(json_files)} JSON, {len(yaml_files)} YAML). All syntax valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
