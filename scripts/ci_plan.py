#!/usr/bin/env python3
"""scripts/ci_plan.py

Deterministic CI Execution Planner and Component/Requirement Impact Analyser.
Models the repository dependency graph across C++ (#include and header-source pairs),
pybind11 symbol bindings (bindings.cpp -> Python AST), and Enhancement Proposal
requirements (Two-Hop Test Impact Analysis + Requirements-Based Testing).

Zero external dependencies (Python standard library only).
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import json
import os
import re
import subprocess
import sys
import time
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import verify_traceability as vt

CPP_INCLUDE_PATTERN = re.compile(r'^\s*#include\s+"([^"]+)"')
GTEST_SUITE_PATTERN = re.compile(
    r"^\s*(?:TEST|TEST_F|TEST_P|TYPED_TEST|TYPED_TEST_P)\s*\(\s*([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*\)"
)
PYBIND_CLASS_PATTERN = re.compile(
    r'py::class_<\s*([A-Za-z0-9_:]+)(?:\s*,[^>]*)?>\s*\(\s*\w+\s*,\s*"([A-Za-z0-9_]+)"'
)
PYBIND_DEF_PATTERN = re.compile(
    r'\.def\(\s*"([A-Za-z0-9_]+)"\s*,\s*&([A-Za-z0-9_]+)::'
)
PYBIND_ATTR_ALIAS_PATTERN = re.compile(
    r'\.attr\(\s*"([A-Za-z0-9_]+)"\s*\)\s*=\s*\w+\.attr\(\s*"([A-Za-z0-9_]+)"\s*\)'
)


@dataclasses.dataclass
class BuildsPlan:
    """Boolean execution flags for downstream CI jobs."""

    validate_docs: bool = False
    lint_cpp: bool = False
    test_ubuntu: bool = False
    test_macos: bool = False
    test_windows: bool = False

    def to_dict(self) -> dict[str, bool]:
        return dataclasses.asdict(self)


@dataclasses.dataclass
class Plan:
    """Deterministic execution plan consumed by CI runners and the Smart Gatekeeper."""

    base_ref: str
    has_changes: bool
    is_full_run: bool
    changed_files: list[str]
    affected_components: list[str]
    affected_cpp_test_files: list[str]
    affected_cpp_suites: list[str]
    ctest_filter: str
    affected_python_test_files: list[str]
    needs_python_bindings_build: bool
    impacted_requirements: list[str]
    claimed_requirements: list[str]
    requirement_triggers: dict[str, list[str]]
    builds: BuildsPlan
    required_checks: list[str]

    def to_dict(self) -> dict[str, Any]:
        data = dataclasses.asdict(self)
        data["builds"] = self.builds.to_dict()
        return data


@dataclasses.dataclass
class DependencyGraph:
    """In-process reverse dependency graph across C++ and Python files."""

    # Map from target file (repo-relative POSIX) -> set of files that directly depend on it
    reverse_deps: dict[str, set[str]]
    # Map from C++ test file -> set of GoogleTest suite names defined in that file
    cpp_test_suites: dict[str, set[str]]
    # All discovered C++ test files
    all_cpp_test_files: set[str]
    # All discovered Python test files
    all_python_test_files: set[str]
    # Python test files that import the compiled `rfl` extension
    python_rfl_test_files: set[str]


def ref_exists(ref: str, repo_root: Path) -> bool:
    """Return True if the given Git reference exists in the repository."""
    if not ref:
        return False
    res = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", ref],
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    return res.returncode == 0


def resolve_base_ref(base_ref_flag: str, repo_root: Path) -> str:
    """Resolve the target Git base reference from CLI flags or CI environment variables."""
    candidates: list[str] = []
    if base_ref_flag:
        candidates.append(base_ref_flag)

    env_base = os.environ.get("CI_BASE_REF") or os.environ.get("GITHUB_BASE_REF")
    if env_base:
        if not env_base.startswith("origin/"):
            candidates.append(f"origin/{env_base}")
        candidates.append(env_base)

    candidates.extend(["origin/main", "main", "HEAD"])
    for candidate in candidates:
        if ref_exists(candidate, repo_root):
            return candidate
    return "HEAD"


def is_full_run(run_all_flag: bool) -> bool:
    """Return True if full-repository regression mode is active."""
    if run_all_flag:
        return True
    if os.environ.get("GITHUB_EVENT_NAME") == "workflow_dispatch":
        return True
    env_val = os.environ.get("CI_RUN_ALL", "").strip().lower()
    return env_val in ("1", "true", "yes")


def _is_ignored_path(rel_posix: str) -> bool:
    """Return True if a path belongs to build outputs or caches that should be ignored."""
    parts = rel_posix.split("/")
    ignored_prefixes = {
        "build",
        "dist",
        "wheelhouse",
        ".venv",
        ".git",
        ".worktrees",
        "node_modules",
    }
    if parts[0] in ignored_prefixes:
        return True
    if "__pycache__" in parts or ".pytest_cache" in parts or ".ruff_cache" in parts:
        return True
    return False


def collect_all_changed_files(
    base_ref: str, repo_root: Path, include_working_tree: bool = True
) -> list[str]:
    """Collect changed files between base_ref and HEAD, plus working tree modifications."""
    files_set: set[str] = set()

    if base_ref and ref_exists(base_ref, repo_root):
        for diff_range in (f"{base_ref}...HEAD", f"{base_ref}..HEAD"):
            res = subprocess.run(
                ["git", "diff", "--name-only", diff_range],
                cwd=repo_root,
                capture_output=True,
                text=True,
            )
            if res.returncode == 0:
                for line in res.stdout.splitlines():
                    clean = line.strip().replace("\\", "/")
                    if clean and not _is_ignored_path(clean):
                        files_set.add(clean)
                break

    if include_working_tree:
        for cmd in (
            ["git", "diff", "--name-only", "HEAD"],
            ["git", "diff", "--name-only"],
            ["git", "diff", "--name-only", "--cached"],
            ["git", "ls-files", "--others", "--exclude-standard"],
        ):
            res = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
            if res.returncode == 0:
                for line in res.stdout.splitlines():
                    clean = line.strip().replace("\\", "/")
                    if clean and not _is_ignored_path(clean):
                        files_set.add(clean)

    return sorted(files_set)


def _resolve_cpp_include(include_str: str, including_file: Path, repo_root: Path) -> str | None:
    """Resolve a C++ #include string to a repository-relative POSIX path if it exists."""
    search_bases = [
        including_file.parent,
        repo_root / "src" / "core",
        repo_root / "src" / "legacy" / "include",
        repo_root / "src" / "legacy" / "source",
        repo_root,
    ]
    for base in search_bases:
        candidate = (base / include_str).resolve()
        if candidate.is_file():
            try:
                return candidate.relative_to(repo_root).as_posix()
            except ValueError:
                continue
    return None


def parse_pybind11_symbol_map(repo_root: Path) -> dict[str, str]:
    """Parse src/python_bindings/bindings.cpp to map Python `rfl.<Symbol>` names to C++ headers.

    This bridges Python tests to the exact C++ headers they exercise without
    requiring manual comment directives.
    """
    bindings_file = repo_root / "src" / "python_bindings" / "bindings.cpp"
    if not bindings_file.exists():
        return {}

    try:
        content = bindings_file.read_text(encoding="utf-8")
    except Exception:
        return {}

    # 1. Map included C++ class/stem names -> repo-relative header path
    stem_to_header: dict[str, str] = {}
    for line in content.splitlines():
        m = CPP_INCLUDE_PATTERN.match(line)
        if m:
            resolved = _resolve_cpp_include(m.group(1), bindings_file, repo_root)
            if resolved:
                stem = Path(resolved).stem
                stem_to_header[stem] = resolved

    symbol_to_header: dict[str, str] = {}

    # 2. Match py::class_<CppClass, ...>(m, "PythonName")
    for m in PYBIND_CLASS_PATTERN.finditer(content):
        cpp_type = m.group(1).split("::")[-1]
        py_name = m.group(2)
        if cpp_type in stem_to_header:
            symbol_to_header[py_name] = stem_to_header[cpp_type]

    # 3. Match m.def("python_func", &CppClass::method)
    for m in PYBIND_DEF_PATTERN.finditer(content):
        py_func = m.group(1)
        cpp_scope = m.group(2).split("::")[-1]
        if cpp_scope in stem_to_header:
            symbol_to_header[py_func] = stem_to_header[cpp_scope]

    # 4. Match module attribute aliases: m.attr("GslRng") = m.attr("StdRng")
    for m in PYBIND_ATTR_ALIAS_PATTERN.finditer(content):
        alias_name = m.group(1)
        target_name = m.group(2)
        if target_name in symbol_to_header:
            symbol_to_header[alias_name] = symbol_to_header[target_name]

    return symbol_to_header


def _analyse_python_ast_dependencies(
    py_file: Path,
    repo_root: Path,
    pybind_symbol_map: dict[str, str],
) -> tuple[set[str], bool]:
    """Inspect a Python file's AST to determine its local module and C++ pybind11 dependencies.

    Returns (set_of_repo_relative_dependencies, imports_rfl_extension).
    """
    deps: set[str] = set()
    imports_rfl = False

    try:
        content = py_file.read_text(encoding="utf-8")
        tree = ast.parse(content, filename=str(py_file))
    except Exception:
        return deps, False

    rfl_aliases: set[str] = set()
    used_rfl_symbols: set[str] = set()
    imported_module_names: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                top_mod = alias.name.split(".")[0]
                imported_module_names.add(top_mod)
                if top_mod == "rfl":
                    imports_rfl = True
                    rfl_aliases.add(alias.asname or "rfl")  # cspell:ignore asname
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                top_mod = node.module.split(".")[0]
                imported_module_names.add(top_mod)
                if top_mod == "rfl":
                    imports_rfl = True
                    for alias in node.names:
                        used_rfl_symbols.add(alias.name)
        elif isinstance(node, ast.Attribute):
            if isinstance(node.value, ast.Name) and node.value.id in rfl_aliases:
                used_rfl_symbols.add(node.attr)

    # Resolve imports from scripts/ or sibling modules
    search_dirs = [py_file.parent, repo_root / "scripts"]
    for mod_name in imported_module_names:
        if mod_name == "rfl":
            continue
        for search_dir in search_dirs:
            candidate = (search_dir / f"{mod_name}.py").resolve()
            if candidate.is_file():
                try:
                    deps.add(candidate.relative_to(repo_root).as_posix())
                except ValueError:
                    pass

    if imports_rfl:
        init_py = repo_root / "src" / "python_bindings" / "rfl" / "__init__.py"
        if init_py.exists():
            deps.add("src/python_bindings/rfl/__init__.py")

        matched_specific_headers: set[str] = set()
        for sym in used_rfl_symbols:
            if sym in pybind_symbol_map:
                matched_specific_headers.add(pybind_symbol_map[sym])

        if matched_specific_headers:
            deps.update(matched_specific_headers)
        else:
            # Generic `import rfl` reflection test (e.g. test_api_surface.py) depends on bindings.cpp
            deps.add("src/python_bindings/bindings.cpp")

    return deps, imports_rfl


def build_dependency_graph(repo_root: Path) -> DependencyGraph:
    """Build the complete C++ (#include + header/source) and Python (AST + pybind11) graph."""
    reverse_deps: dict[str, set[str]] = {}
    cpp_test_suites: dict[str, set[str]] = {}
    all_cpp_test_files: set[str] = set()
    all_python_test_files: set[str] = set()
    python_rfl_test_files: set[str] = set()

    def add_edge(dependency: str, dependent: str) -> None:
        if dependency != dependent:
            reverse_deps.setdefault(dependency, set()).add(dependent)

    # 1. Scan C++ files across src/, examples/, playground/, tests/
    cpp_dirs = [
        repo_root / "src",
        repo_root / "examples",
        repo_root / "playground",
        repo_root / "tests",
    ]
    for cpp_dir in cpp_dirs:
        if not cpp_dir.exists():
            continue
        for ext in ("*.cpp", "*.hpp", "*.h", "*.cc"):
            for cpp_file in cpp_dir.rglob(ext):
                if "OldCode" in cpp_file.parts or "performance" in cpp_file.parts:
                    continue
                rel_path = cpp_file.relative_to(repo_root).as_posix()
                is_test = vt.is_test_file_path(rel_path)

                try:
                    lines = cpp_file.read_text(encoding="utf-8").splitlines()
                except Exception:
                    continue

                suites_in_file: set[str] = set()
                for line in lines:
                    inc_match = CPP_INCLUDE_PATTERN.match(line)
                    if inc_match:
                        inc_target = _resolve_cpp_include(inc_match.group(1), cpp_file, repo_root)
                        if inc_target:
                            # If inc_target changes, rel_path is affected
                            add_edge(inc_target, rel_path)
                            # Header-implementation pairing: if Foo.cpp includes Foo.hpp (and is not a test),
                            # modifying Foo.cpp affects Foo.hpp and all consumers of Foo.hpp.
                            if (
                                not is_test
                                and cpp_file.suffix in (".cpp", ".cc")
                                and Path(inc_target).stem == cpp_file.stem
                            ):
                                add_edge(rel_path, inc_target)

                    if is_test:
                        gm = GTEST_SUITE_PATTERN.match(line)
                        if gm:
                            suites_in_file.add(gm.group(1))

                if is_test and suites_in_file:
                    all_cpp_test_files.add(rel_path)
                    cpp_test_suites[rel_path] = suites_in_file

    # 2. Parse pybind11 symbol bindings in src/python_bindings/bindings.cpp
    pybind_symbol_map = parse_pybind11_symbol_map(repo_root)

    # 3. Scan Python files across src/python_bindings/, scripts/, tests/
    py_dirs = [
        repo_root / "src" / "python_bindings",
        repo_root / "scripts",
        repo_root / "tests",
    ]
    for py_dir in py_dirs:
        if not py_dir.exists():
            continue
        for py_file in py_dir.rglob("*.py"):
            rel_path = py_file.relative_to(repo_root).as_posix()
            if _is_ignored_path(rel_path):
                continue

            is_test = vt.is_test_file_path(rel_path) and py_file.name.startswith("test_")
            if is_test:
                all_python_test_files.add(rel_path)

            deps, imports_rfl = _analyse_python_ast_dependencies(
                py_file, repo_root, pybind_symbol_map
            )
            for dep in deps:
                add_edge(dep, rel_path)

            if is_test and imports_rfl:
                python_rfl_test_files.add(rel_path)

    return DependencyGraph(
        reverse_deps=reverse_deps,
        cpp_test_suites=cpp_test_suites,
        all_cpp_test_files=all_cpp_test_files,
        all_python_test_files=all_python_test_files,
        python_rfl_test_files=python_rfl_test_files,
    )


def extract_modified_ep_requirements(
    base_ref: str,
    repo_root: Path,
    req_map: dict[str, vt.Requirement],
    include_working_tree: bool = True,
) -> set[str]:
    """Identify requirement IDs whose table rows in docs/eps/ep-*.md were added or edited."""
    modified_requirements: set[str] = set()
    diff_commands: list[list[str]] = []

    if base_ref and ref_exists(base_ref, repo_root):
        diff_commands.append(["git", "diff", "-U0", f"{base_ref}...HEAD", "--", "docs/eps/"])
    if include_working_tree:
        diff_commands.append(["git", "diff", "-U0", "HEAD", "--", "docs/eps/"])

    current_ep_file = ""
    for cmd in diff_commands:
        res = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
        if res.returncode != 0:
            continue
        for line in res.stdout.splitlines():
            if line.startswith("+++ b/"):
                current_ep_file = line[6:].strip()
                continue
            if line.startswith(("+++", "---", "@@")):
                continue
            if line.startswith(("+", "-")):
                # Only trigger on table rows defining a requirement (first cell has REQ-...)
                stripped = line[1:].strip()
                if not stripped.startswith("|"):
                    continue
                cells = [c.strip() for c in stripped.strip("|").split("|")]
                if not cells:
                    continue
                first_cell = cells[0]
                for raw_id in vt.REQ_TAG_SEARCH_PATTERN.findall(first_cell):
                    can_id = vt.canonicalize_req(raw_id, current_ep_file)
                    if can_id in req_map:
                        modified_requirements.add(can_id)

    return modified_requirements


def extract_claimed_trailer_requirements(base_ref: str, repo_root: Path) -> set[str]:
    """Extract requirement IDs claimed via 'Implements: REQ-...' in commits or PR description."""
    texts_to_scan: list[str] = []
    if base_ref and ref_exists(base_ref, repo_root):
        res = subprocess.run(
            ["git", "log", f"{base_ref}..HEAD", "--pretty=format:%B---COMMIT---"],
            cwd=repo_root,
            capture_output=True,
            text=True,
        )
        if res.returncode == 0 and res.stdout.strip():
            texts_to_scan.extend(res.stdout.split("---COMMIT---"))

    if not texts_to_scan:
        res = subprocess.run(
            ["git", "log", "-1", "--pretty=format:%B"],
            cwd=repo_root,
            capture_output=True,
            text=True,
        )
        if res.returncode == 0:
            texts_to_scan.append(res.stdout)

    event_path = os.environ.get("GITHUB_EVENT_PATH")
    if event_path and Path(event_path).exists():
        try:
            event_data = json.loads(Path(event_path).read_text(encoding="utf-8"))
            pr_body = event_data.get("pull_request", {}).get("body")
            if pr_body:
                texts_to_scan.append(pr_body)
        except Exception:
            pass

    claimed: set[str] = set()
    trailer_regex = re.compile(r"(?im)^[\s*-]*Implements:\s*(.+)$")
    for txt in texts_to_scan:
        for m in trailer_regex.finditer(txt):
            for token in vt.REQ_TAG_SEARCH_PATTERN.findall(m.group(1)):
                claimed.add(vt.canonicalize_req(token))

    return claimed


def compute_plan_from_files(
    changed_files: list[str],
    base_ref: str,
    repo_root: Path,
    eps_dir: Path | None = None,
    full_run: bool = False,
    ep_diff_requirements: set[str] | None = None,
    claimed_requirements: set[str] | None = None,
) -> Plan:
    """Compute the deterministic CI execution plan from a list of changed files."""
    resolved_eps_dir = (eps_dir or (repo_root / "docs" / "eps")).resolve()
    graph = build_dependency_graph(repo_root)
    req_map, _ = vt.collect_all_requirements(resolved_eps_dir, repo_root)
    all_test_mappings = vt.collect_all_tests(repo_root)

    if full_run:
        all_cpp_files = sorted(graph.all_cpp_test_files)
        all_suites = sorted({s for suites in graph.cpp_test_suites.values() for s in suites})
        all_py_files = sorted(graph.all_python_test_files)
        implemented_requirements = sorted(
            r.req_id for r in req_map.values() if r.is_implemented and r.is_unit_testable
        )
        builds = BuildsPlan(
            validate_docs=True,
            lint_cpp=True,
            test_ubuntu=True,
            test_macos=True,
            test_windows=True,
        )
        return Plan(
            base_ref=base_ref,
            has_changes=True,
            is_full_run=True,
            changed_files=changed_files,
            affected_components=["all"],
            affected_cpp_test_files=all_cpp_files,
            affected_cpp_suites=all_suites,
            ctest_filter="",
            affected_python_test_files=all_py_files,
            needs_python_bindings_build=bool(graph.python_rfl_test_files),
            impacted_requirements=implemented_requirements,
            claimed_requirements=sorted(claimed_requirements or set()),
            requirement_triggers={r: ["full-run"] for r in implemented_requirements},
            builds=builds,
            required_checks=[
                "validate_docs",
                "lint",
                "test_ubuntu",
                "test_macos",
                "test_windows",
            ],
        )

    if not changed_files and not claimed_requirements and not ep_diff_requirements:
        return Plan(
            base_ref=base_ref,
            has_changes=False,
            is_full_run=False,
            changed_files=[],
            affected_components=[],
            affected_cpp_test_files=[],
            affected_cpp_suites=[],
            ctest_filter="",
            affected_python_test_files=[],
            needs_python_bindings_build=False,
            impacted_requirements=[],
            claimed_requirements=[],
            requirement_triggers={},
            builds=BuildsPlan(),
            required_checks=[],
        )

    # Detect build-wide or workflow-wide changes
    cpp_build_wide = any(
        f in ("CMakeLists.txt", "src/core/CMakeLists.txt", "src/core/tests/CMakeLists.txt")
        or f.startswith("cmake/")
        or f.startswith(".github/actions/")
        for f in changed_files
    )
    py_build_wide = any(
        f in ("pyproject.toml", "src/python_bindings/CMakeLists.txt") for f in changed_files
    )
    workflow_changed = any(
        f.startswith(".github/workflows/") or f == "mise.toml" for f in changed_files
    )
    cpp_source_or_header_changed = any(
        f.endswith((".cpp", ".hpp", ".h", ".cc")) and not vt.is_test_file_path(f)
        for f in changed_files
    )
    any_cpp_file_changed = any(
        f.endswith((".cpp", ".hpp", ".h", ".cc")) or f == ".clang-format" for f in changed_files
    )

    # --- HOP 1: Reverse Transitive Dependency Closure (Code -> Affected Files & Tests) ---
    affected_files: set[str] = set(changed_files)
    queue = list(changed_files)

    # If bindings.cpp itself was directly modified, all Python `rfl` tests are affected
    if "src/python_bindings/bindings.cpp" in affected_files:
        for py_rfl_test in graph.python_rfl_test_files:
            if py_rfl_test not in affected_files:
                affected_files.add(py_rfl_test)
                queue.append(py_rfl_test)

    # If pyproject.toml or src/python_bindings/CMakeLists.txt changed, all Python tests are affected
    if py_build_wide:
        for py_test in graph.all_python_test_files:
            if py_test not in affected_files:
                affected_files.add(py_test)
                queue.append(py_test)

    if cpp_build_wide:
        for cpp_test in graph.all_cpp_test_files:
            if cpp_test not in affected_files:
                affected_files.add(cpp_test)
                queue.append(cpp_test)

    while queue:
        curr = queue.pop(0)
        for dependent in graph.reverse_deps.get(curr, set()):
            if dependent not in affected_files:
                affected_files.add(dependent)
                queue.append(dependent)

    affected_cpp_test_files: set[str] = {
        f for f in affected_files if f in graph.all_cpp_test_files
    }
    affected_python_test_files: set[str] = {
        f for f in affected_files if f in graph.all_python_test_files
    }

    # --- HOP 2: Impacted Requirements Discovery (Tests + EP Diffs + Trailers -> REQs) ---
    req_triggers: dict[str, set[str]] = {}

    # 2a. From affected component tests (Two-Hop inference without source code tags)
    for mapping in all_test_mappings:
        if (
            mapping.source_file in affected_cpp_test_files
            or mapping.source_file in affected_python_test_files
        ):
            if mapping.req_id in req_map:
                req = req_map[mapping.req_id]
                if req.is_implemented and req.is_unit_testable:
                    req_triggers.setdefault(mapping.req_id, set()).add("component-test")

    # 2b. From EP specification table diffs in docs/eps/ep-*.md
    if ep_diff_requirements:
        for req_id in ep_diff_requirements:
            if req_id in req_map:
                req = req_map[req_id]
                if req.is_implemented and req.is_unit_testable:
                    req_triggers.setdefault(req_id, set()).add("ep-spec-diff")

    # 2c. From explicit Git trailers ("Implements: REQ-...")
    if claimed_requirements:
        for req_id in claimed_requirements:
            req_triggers.setdefault(req_id, set()).add("git-trailer")

    impacted_requirements = sorted(req_triggers.keys())

    # --- HOP 3: Requirement Test Closure (Impacted REQs -> All Mapped Tests) ---
    for mapping in all_test_mappings:
        if mapping.req_id in req_triggers:
            if mapping.source_file in graph.all_cpp_test_files:
                affected_cpp_test_files.add(mapping.source_file)
            elif mapping.source_file in graph.all_python_test_files:
                affected_python_test_files.add(mapping.source_file)

    # Collect C++ GoogleTest suite names for all affected C++ test files
    affected_cpp_suites: set[str] = set()
    for tf in affected_cpp_test_files:
        affected_cpp_suites.update(graph.cpp_test_suites.get(tf, set()))

    # If workflow files changed and no specific C++ test suites were otherwise affected,
    # run all C++ suites to validate the workflow change across platforms.
    if cpp_build_wide or (workflow_changed and not affected_cpp_suites):
        affected_cpp_test_files = set(graph.all_cpp_test_files)
        affected_cpp_suites = {s for suites in graph.cpp_test_suites.values() for s in suites}
        ctest_filter = ""
    elif affected_cpp_suites:
        sorted_suites = sorted(affected_cpp_suites)
        ctest_filter = f"^({'|'.join(sorted_suites)})\\."
    else:
        ctest_filter = ""

    needs_python_bindings_build = any(
        tf in graph.python_rfl_test_files for tf in affected_python_test_files
    )

    # Determine non-test affected components for summary reporting
    affected_components = sorted(
        f
        for f in affected_files
        if not vt.is_test_file_path(f)
        and (
            f.startswith(("src/", "scripts/", "examples/", "playground/", "cmake/"))
            or f in ("CMakeLists.txt", "pyproject.toml", "mise.toml")
        )
    )

    # Determine which CI jobs need to execute
    run_validate_docs = True
    run_lint_cpp = any_cpp_file_changed or cpp_build_wide or workflow_changed
    run_test_ubuntu = (
        bool(affected_cpp_suites)
        or cpp_source_or_header_changed
        or cpp_build_wide
        or workflow_changed
    )
    # Run macOS and Windows when C++ production code, CMake, or workflows change
    run_cross_platform = cpp_source_or_header_changed or cpp_build_wide or workflow_changed
    run_test_macos = run_cross_platform
    run_test_windows = run_cross_platform

    builds = BuildsPlan(
        validate_docs=run_validate_docs,
        lint_cpp=run_lint_cpp,
        test_ubuntu=run_test_ubuntu,
        test_macos=run_test_macos,
        test_windows=run_test_windows,
    )

    required_checks: list[str] = []
    if builds.validate_docs:
        required_checks.append("validate_docs")
    if builds.lint_cpp:
        required_checks.append("lint")
    if builds.test_ubuntu:
        required_checks.append("test_ubuntu")
    if builds.test_macos:
        required_checks.append("test_macos")
    if builds.test_windows:
        required_checks.append("test_windows")

    return Plan(
        base_ref=base_ref,
        has_changes=True,
        is_full_run=False,
        changed_files=sorted(changed_files),
        affected_components=affected_components,
        affected_cpp_test_files=sorted(affected_cpp_test_files),
        affected_cpp_suites=sorted(affected_cpp_suites),
        ctest_filter=ctest_filter,
        affected_python_test_files=sorted(affected_python_test_files),
        needs_python_bindings_build=needs_python_bindings_build,
        impacted_requirements=impacted_requirements,
        claimed_requirements=sorted(claimed_requirements or set()),
        requirement_triggers={k: sorted(v) for k, v in sorted(req_triggers.items())},
        builds=builds,
        required_checks=required_checks,
    )


def compute_plan(
    base_ref_flag: str = "",
    repo_root: Path | None = None,
    eps_dir: Path | None = None,
    run_all: bool = False,
    include_working_tree: bool = True,
) -> Plan:
    """Resolve Git state and compute the full CI and requirement validation plan."""
    root = (repo_root or SCRIPT_DIR.parent).resolve()
    resolved_eps_dir = (eps_dir or (root / "docs" / "eps")).resolve()
    base_ref = resolve_base_ref(base_ref_flag, root)
    full = is_full_run(run_all)

    changed_files = collect_all_changed_files(
        base_ref, root, include_working_tree=include_working_tree
    )
    req_map, _ = vt.collect_all_requirements(resolved_eps_dir, root)
    ep_diff_requirements = extract_modified_ep_requirements(
        base_ref, root, req_map, include_working_tree=include_working_tree
    )
    claimed_requirements = extract_claimed_trailer_requirements(base_ref, root)

    return compute_plan_from_files(
        changed_files=changed_files,
        base_ref=base_ref,
        repo_root=root,
        eps_dir=resolved_eps_dir,
        full_run=full,
        ep_diff_requirements=ep_diff_requirements,
        claimed_requirements=claimed_requirements,
    )


def run_unittest_with_junit_xml(test_files: list[str], repo_root: Path, xml_path: Path | None) -> int:
    """Execute standard-library Python unittest files and write a JUnit XML report."""
    if not test_files:
        return 0

    suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    for rel_file in test_files:
        abs_file = (repo_root / rel_file).resolve()
        if not abs_file.exists():
            continue
        parent_dir = str(abs_file.parent)
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)
        discovered = loader.discover(start_dir=parent_dir, pattern=abs_file.name)
        suite.addTests(discovered)

    class XMLTrackingResult(unittest.TextTestResult):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.records: list[tuple[unittest.TestCase, str, str, float]] = []
            self._start_time: float = 0.0

        def startTest(self, test: unittest.TestCase) -> None:
            self._start_time = time.perf_counter()
            super().startTest(test)

        def addSuccess(self, test: unittest.TestCase) -> None:
            elapsed = time.perf_counter() - self._start_time
            self.records.append((test, "PASSED", "", elapsed))
            super().addSuccess(test)

        def addFailure(self, test: unittest.TestCase, err: Any) -> None:
            elapsed = time.perf_counter() - self._start_time
            msg = self._exc_info_to_string(err, test)
            self.records.append((test, "FAILED", msg, elapsed))
            super().addFailure(test, err)

        def addError(self, test: unittest.TestCase, err: Any) -> None:
            elapsed = time.perf_counter() - self._start_time
            msg = self._exc_info_to_string(err, test)
            self.records.append((test, "ERROR", msg, elapsed))
            super().addError(test, err)

        def addSkip(self, test: unittest.TestCase, reason: str) -> None:
            elapsed = time.perf_counter() - self._start_time
            self.records.append((test, "SKIPPED", reason, elapsed))
            super().addSkip(test, reason)

    runner = unittest.TextTestRunner(verbosity=2, resultclass=XMLTrackingResult)
    result: XMLTrackingResult = runner.run(suite)  # type: ignore[assignment]

    if xml_path is not None:
        xml_path.parent.mkdir(parents=True, exist_ok=True)
        ts_el = ET.Element(
            "testsuite",
            attrib={
                "name": "python_unittest",
                "tests": str(len(result.records)),
                "failures": str(len(result.failures)),
                "errors": str(len(result.errors)),
                "skipped": str(len(result.skipped)),
            },
        )
        for tc, status, message, elapsed in result.records:
            tc_id = tc.id()
            if "." in tc_id:
                classname, name = tc_id.rsplit(".", 1)
            else:
                classname, name = "unittest", tc_id

            tc_el = ET.SubElement(
                ts_el,
                "testcase",
                attrib={
                    "classname": classname,
                    "name": name,
                    "time": f"{elapsed:.4f}",
                },
            )
            if status == "FAILED":
                fail_el = ET.SubElement(tc_el, "failure", attrib={"message": "AssertionError"})
                fail_el.text = message
            elif status == "ERROR":
                err_el = ET.SubElement(tc_el, "error", attrib={"message": "ExecutionError"})
                err_el.text = message
            elif status == "SKIPPED":
                ET.SubElement(tc_el, "skipped", attrib={"message": message})

        tree = ET.ElementTree(ts_el)
        ET.indent(tree, space="  ")
        tree.write(xml_path, encoding="utf-8", xml_declaration=True)

    return 0 if result.wasSuccessful() else 1


def print_human_summary(plan: Plan) -> None:
    """Print a mobile-friendly human summary of the execution plan."""
    print(f"==> RFL CI Execution Plan (Base: {plan.base_ref}, FullRun: {plan.is_full_run})")
    print(f"  • Has Changes:            {plan.has_changes}")
    print(f"  • Changed Files:          {len(plan.changed_files)}")
    print(f"  • Affected Components:    {len(plan.affected_components)}")
    for comp in plan.affected_components:
        print(f"      - {comp}")

    print(f"  • Affected C++ Suites:    {len(plan.affected_cpp_suites)}")
    for suite in plan.affected_cpp_suites:
        print(f"      - {suite}")
    if plan.ctest_filter:
        print(f"  • CTest Filter:           {plan.ctest_filter}")

    print(f"  • Affected Python Tests:  {len(plan.affected_python_test_files)}")
    for py_tf in plan.affected_python_test_files:
        print(f"      - {py_tf}")

    print(f"  • Impacted Requirements:  {len(plan.impacted_requirements)}")
    for req_id in plan.impacted_requirements:
        reasons = ", ".join(plan.requirement_triggers.get(req_id, []))
        print(f"      - {req_id} ({reasons})")

    print(f"  • Scheduled CI Checks:    {', '.join(plan.required_checks) or 'none'}")


def write_step_summary(plan: Plan) -> None:
    """Append the execution plan summary to $GITHUB_STEP_SUMMARY when running in GitHub Actions."""
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    lines: list[str] = [
        "### 🧭 RFL CI Execution Plan",
        "",
        f"- **Base Reference:** `{plan.base_ref}`",
        f"- **Full Regression Mode:** `{plan.is_full_run}`",
        f"- **Scheduled Checks:** {', '.join(f'`{c}`' for c in plan.required_checks) or '*None*'}",
        "",
    ]

    if plan.affected_cpp_suites:
        suites_str = ", ".join(f"`{s}`" for s in plan.affected_cpp_suites)
        lines.append(
            f"- **Affected C++ Test Suites ({len(plan.affected_cpp_suites)}):** {suites_str}"
        )
    else:
        lines.append("- **Affected C++ Test Suites:** *None*")

    if plan.affected_python_test_files:
        py_str = ", ".join(f"`{Path(f).name}`" for f in plan.affected_python_test_files)
        lines.append(
            f"- **Affected Python Test Files ({len(plan.affected_python_test_files)}):** {py_str}"
        )
    else:
        lines.append("- **Affected Python Test Files:** *None*")

    if plan.impacted_requirements:
        req_items = [
            f"`{r}` ({'+'.join(plan.requirement_triggers.get(r, []))})"
            for r in plan.impacted_requirements
        ]
        lines.append(
            f"- **Impacted Requirements ({len(plan.impacted_requirements)}):** {', '.join(req_items)}"
        )
    else:
        lines.append("- **Impacted Requirements:** *None*")

    lines.append("")
    with open(summary_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def write_github_output(plan: Plan) -> None:
    """Export execution plan keys to $GITHUB_OUTPUT when running in GitHub Actions."""
    output_path = os.environ.get("GITHUB_OUTPUT")
    if not output_path:
        return
    lines = [
        f"has_changes={'true' if plan.has_changes else 'false'}",
        f"validate_docs={'true' if plan.builds.validate_docs else 'false'}",
        f"lint_cpp={'true' if plan.builds.lint_cpp else 'false'}",
        f"test_ubuntu={'true' if plan.builds.test_ubuntu else 'false'}",
        f"test_macos={'true' if plan.builds.test_macos else 'false'}",
        f"test_windows={'true' if plan.builds.test_windows else 'false'}",
        f"ctest_filter={plan.ctest_filter}",
        f"required_checks={json.dumps(plan.required_checks)}",
    ]
    with open(output_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def cmd_plan(args: argparse.Namespace, repo_root: Path) -> int:
    """Execute the `plan` sub-command."""
    plan = compute_plan(
        base_ref_flag=args.base_ref,
        repo_root=repo_root,
        run_all=args.all,
        include_working_tree=not bool(os.environ.get("CI")),
    )
    payload = json.dumps(plan.to_dict(), indent=2)

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(payload + "\n", encoding="utf-8")

    write_github_output(plan)
    write_step_summary(plan)

    if args.summary and not args.json:
        print_human_summary(plan)
    else:
        print(payload)
    return 0


def cmd_affected(args: argparse.Namespace, repo_root: Path) -> int:
    """Execute the `affected` sub-command."""
    plan = compute_plan(
        base_ref_flag=args.base_ref,
        repo_root=repo_root,
        run_all=args.all,
        include_working_tree=True,
    )
    if args.json:
        print(json.dumps(plan.to_dict(), indent=2))
    else:
        print_human_summary(plan)
    return 0


def cmd_test_affected(args: argparse.Namespace, repo_root: Path) -> int:
    """Execute the `test-affected` sub-command on affected C++ and Python test suites."""
    plan = compute_plan(
        base_ref_flag=args.base_ref,
        repo_root=repo_root,
        run_all=args.all,
        include_working_tree=not bool(os.environ.get("CI")),
    )
    print_human_summary(plan)

    exit_code = 0

    # 1. Run affected C++ tests via CTest
    if not args.only_python:
        if plan.affected_cpp_suites or plan.is_full_run:
            build_dir = repo_root / args.build_dir
            if not args.no_build:
                if not build_dir.exists():
                    print(f"--> Configuring and building C++ targets in {build_dir}...")
                    subprocess.run(
                        ["cmake", "-B", str(build_dir), "-S", str(repo_root), "-G", "Ninja"],
                        cwd=repo_root,
                        check=True,
                    )
                build_res = subprocess.run(["cmake", "--build", str(build_dir)], cwd=repo_root)
                if build_res.returncode != 0 and not (build_dir / "CTestTestfile.cmake").exists():
                    return build_res.returncode

            ctest_cmd = [
                "ctest",
                "--test-dir",
                str(build_dir),
                "--output-on-failure",
            ]
            if plan.ctest_filter:
                ctest_cmd.extend(["-R", plan.ctest_filter])
            if args.cpp_xml:
                xml_abs = (repo_root / args.cpp_xml).resolve()
                xml_abs.parent.mkdir(parents=True, exist_ok=True)
                ctest_cmd.extend(["--output-junit", str(xml_abs)])

            print(f"--> Running CTest: {' '.join(ctest_cmd)}")
            res = subprocess.run(ctest_cmd, cwd=repo_root)
            if res.returncode != 0:
                exit_code = res.returncode
        else:
            print("--> No C++ test suites affected. Skipping CTest.")

    # 2. Run affected Python tests
    if not args.only_cpp:
        py_files = plan.affected_python_test_files
        if py_files:
            stdlib_tests = [f for f in py_files if f.startswith("scripts/tests/")]
            rfl_tests = [f for f in py_files if not f.startswith("scripts/tests/")]

            py_xml_path = (repo_root / args.py_xml).resolve() if args.py_xml else None

            if stdlib_tests:
                print(f"--> Running affected script unit tests: {', '.join(stdlib_tests)}")
                rc = run_unittest_with_junit_xml(stdlib_tests, repo_root, py_xml_path)
                if rc != 0:
                    exit_code = rc

            if rfl_tests and not os.environ.get("RFL_SKIP_PYBIND_TESTS"):
                # Check if uv or pytest is available to run rfl binding tests
                uv_bin = subprocess.run(["which", "uv"], capture_output=True, text=True)
                if uv_bin.returncode == 0 and not os.environ.get("CI"):
                    pytest_cmd = ["uv", "run", "--group", "test", "pytest", *rfl_tests]
                    print(f"--> Running affected Python binding tests: {' '.join(pytest_cmd)}")
                    res = subprocess.run(pytest_cmd, cwd=repo_root)
                    if res.returncode != 0:
                        exit_code = res.returncode
        else:
            print("--> No Python test files affected. Skipping Python tests.")

    return exit_code


def cmd_check_milestone() -> int:
    """Verify that a GitHub pull request has an assigned milestone."""
    event_name = os.environ.get("GITHUB_EVENT_NAME", "pull_request")
    if event_name != "pull_request":
        print(f"Skipping milestone check for event: {event_name}")
        return 0

    milestone = ""
    pr_number: str | int = ""
    repo = os.environ.get("GITHUB_REPOSITORY", "")

    event_path = os.environ.get("GITHUB_EVENT_PATH")
    if event_path and Path(event_path).exists():
        try:
            event_data = json.loads(Path(event_path).read_text(encoding="utf-8"))
            pr_data = event_data.get("pull_request") or {}
            pr_number = pr_data.get("number", "")
            milestone_obj = pr_data.get("milestone") or {}
            milestone = (milestone_obj.get("title") or "").strip()
        except Exception:
            pass

    if not milestone and pr_number and repo:
        res = subprocess.run(
            [
                "gh",
                "pr",
                "view",
                str(pr_number),
                "--repo",
                repo,
                "--json",
                "milestone",
                "--jq",
                ".milestone.title // empty",
            ],
            capture_output=True,
            text=True,
        )
        if res.returncode == 0:
            milestone = res.stdout.strip()

    if not milestone:
        pr_label = f"#{pr_number}" if pr_number else "current PR"
        print(
            f"::error::Pull Request {pr_label} does not have an assigned milestone. "
            "All PRs targeting main must belong to an active milestone.",
            file=sys.stderr,
        )
        return 1

    print(f"✅ Pull request is assigned to milestone: {milestone}")
    return 0


def cmd_gate() -> int:
    """Evaluate upstream CI job statuses from $NEEDS_JSON against required_checks."""
    needs_raw = os.environ.get("NEEDS_JSON", "").strip()
    if not needs_raw:
        print("❌ NEEDS_JSON environment variable is empty or unset.", file=sys.stderr)
        return 1

    try:
        needs = json.loads(needs_raw)
    except Exception as exc:
        print(f"❌ Failed to parse NEEDS_JSON: {exc}", file=sys.stderr)
        return 1

    if needs.get("plan", {}).get("result") != "success":
        print("❌ CI Execution Planner failed or was cancelled.", file=sys.stderr)
        return 1

    if needs.get("check_milestone", {}).get("result") in ("failure", "cancelled"):
        print(
            "❌ Milestone verification failed: pull request must have an assigned milestone.",
            file=sys.stderr,
        )
        return 1

    required_raw = needs.get("plan", {}).get("outputs", {}).get("required_checks", "[]")
    try:
        required: list[str] = json.loads(required_raw or "[]")
    except Exception:
        required = []

    failed = [j for j in required if needs.get(j, {}).get("result") != "success"]
    for matrix_job in ("matrix_linux", "matrix_macos", "matrix_windows"):
        if needs.get(matrix_job, {}).get("result") in ("failure", "cancelled"):
            failed.append(matrix_job)

    if failed:
        print(f"❌ Failed or cancelled required checks: {failed}", file=sys.stderr)
        return 1

    print(f"✅ All required checks passed: {required}")
    return 0


def cmd_ctest(args: argparse.Namespace, repo_root: Path) -> int:
    """Execute CTest with optional regex filter, parallel workers, and JUnit XML."""
    build_dir = (repo_root / args.build_dir).resolve()
    if not build_dir.exists():
        print(f"❌ Build directory not found: {build_dir}", file=sys.stderr)
        return 1

    parallel = (
        args.parallel
        or os.environ.get("CTEST_PARALLEL_LEVEL")
        or str(os.cpu_count() or 2)
    )
    build_type = args.build_type or os.environ.get("BUILD_TYPE") or "Release"
    filter_regex = (args.filter or os.environ.get("CTEST_FILTER") or "").strip()
    junit_xml = (args.junit or os.environ.get("CTEST_JUNIT") or "").strip()

    ctest_cmd = [
        "ctest",
        "--test-dir",
        str(build_dir),
        "-C",
        str(build_type),
        "--parallel",
        str(parallel),
        "--output-on-failure",
    ]

    if filter_regex:
        ctest_cmd.extend(["-R", filter_regex])

    if junit_xml:
        xml_path = (repo_root / junit_xml).resolve()
        xml_path.parent.mkdir(parents=True, exist_ok=True)
        ctest_cmd.extend(["--output-junit", str(xml_path)])

    extra = getattr(args, "extra_ctest_args", [])
    if extra:
        if extra[0] == "--":
            extra = extra[1:]
        ctest_cmd.extend(extra)

    print(f"--> Running CTest: {' '.join(ctest_cmd)}")
    res = subprocess.run(ctest_cmd, cwd=repo_root)
    return res.returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="RFL Deterministic CI Execution Planner & Component/Requirement Impact Analyser"
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=SCRIPT_DIR.parent,
        help="Repository root directory",
    )

    sub_parsers = parser.add_subparsers(dest="command")

    # `plan` sub-command
    plan_p = sub_parsers.add_parser(
        "plan", help="Compute deterministic execution plan as JSON for CI runners and gatekeepers"
    )
    plan_p.add_argument("--base-ref", type=str, default="", help="Target Git base ref")
    plan_p.add_argument(
        "-a", "--all", action="store_true", help="Generate full regression plan for all suites"
    )
    plan_p.add_argument(
        "--json", action="store_true", default=False, help="Output plan as formatted JSON"
    )
    plan_p.add_argument(
        "--summary", action="store_true", default=False, help="Print human-readable summary"
    )
    plan_p.add_argument("-o", "--output", type=str, default="", help="Write JSON plan to file")

    # `affected` sub-command
    aff_p = sub_parsers.add_parser(
        "affected",
        help="List affected components, test suites, and requirements (including uncommitted work)",
    )
    aff_p.add_argument("--base-ref", type=str, default="", help="Target Git base ref")
    aff_p.add_argument("-a", "--all", action="store_true", help="List all repository suites")
    aff_p.add_argument("--json", action="store_true", help="Output as JSON")

    # `test-affected` sub-command
    test_p = sub_parsers.add_parser(
        "test-affected",
        help="Run unit tests strictly on affected C++ suites and Python test files",
    )
    test_p.add_argument("--base-ref", type=str, default="", help="Target Git base ref")
    test_p.add_argument("-a", "--all", action="store_true", help="Run all tests")
    test_p.add_argument(
        "--no-build", action="store_true", help="Skip cmake --build step before running CTest"
    )
    test_p.add_argument(
        "--only-python", action="store_true", help="Run only affected Python test files"
    )
    test_p.add_argument(
        "--only-cpp", action="store_true", help="Run only affected C++ test suites"
    )
    test_p.add_argument(
        "--build-dir", type=str, default="build", help="CMake build directory (default: build)"
    )
    test_p.add_argument(
        "--cpp-xml", type=str, default="", help="Output path for CTest JUnit XML report"
    )
    test_p.add_argument(
        "--py-xml", type=str, default="", help="Output path for Python JUnit XML report"
    )

    # `check-milestone` sub-command
    sub_parsers.add_parser(
        "check-milestone", help="Verify that the current GitHub pull request has an assigned milestone"
    )

    # `gate` sub-command
    sub_parsers.add_parser(
        "gate", help="Evaluate upstream CI job statuses from $NEEDS_JSON against required_checks"
    )

    # `ctest` sub-command
    ctest_p = sub_parsers.add_parser(
        "ctest",
        help="Execute CTest with optional regex filter, parallel workers, and JUnit XML",
    )
    ctest_p.add_argument(
        "-R", "--filter", type=str, default="", help="Regular expression filter for test names"
    )
    ctest_p.add_argument(
        "-C", "--build-type", type=str, default="", help="CMake build configuration (Release, Debug)"
    )
    ctest_p.add_argument(
        "-j", "--parallel", type=str, default="", help="Number of parallel test jobs"
    )
    ctest_p.add_argument(
        "--build-dir", type=str, default="build", help="CMake build directory (default: build)"
    )
    ctest_p.add_argument(
        "--junit", type=str, default="", help="Output path for CTest JUnit XML report"
    )
    ctest_p.add_argument(
        "extra_ctest_args",
        nargs=argparse.REMAINDER,
        help="Additional arguments forwarded directly to CTest",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()

    if args.command == "affected":
        return cmd_affected(args, repo_root)
    if args.command == "test-affected":
        return cmd_test_affected(args, repo_root)
    if args.command == "plan":
        return cmd_plan(args, repo_root)
    if args.command == "check-milestone":
        return cmd_check_milestone()
    if args.command == "gate":
        return cmd_gate()
    if args.command == "ctest":
        return cmd_ctest(args, repo_root)

    # Default to `affected` human summary if invoked without sub-command
    default_args = argparse.Namespace(base_ref="", all=False, json=False)
    return cmd_affected(default_args, repo_root)


if __name__ == "__main__":
    sys.exit(main())
