#!/usr/bin/env python3
"""scripts/verify_traceability.py

Specification-Driven Requirements Traceability and Governance Engine.
Verifies compliance of repository specifications, test suites, and pull requests
against Enhancement Proposal (EP) requirements across GoogleTest, Catch2, and Pytest.

Conforms strictly to EP-5 (REQ-TRC-001 through REQ-TRC-008).
Uses standard Python library only (zero external dependencies).
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import hashlib
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

# Global Requirement ID specification (REQ-TRC-001)
REQ_ID_PATTERN = re.compile(r"^REQ-[A-Z]{2,4}-\d{3}$")
REQ_TAG_SEARCH_PATTERN = re.compile(r"\b(REQ-[A-Za-z0-9_-]+)\b")

# Legacy requirement aliases mapping (prefix_or_none, legacy_id) -> canonical_id
# Resolves collisions (such as REQ-001 in EP-1 vs EP-4) without modifying main test files.
LEGACY_REQ_ALIASES: dict[tuple[str | None, str], str] = {
    # Core architecture disambiguation (EP-1)
    ("src/core/tests", "REQ-001"): "REQ-ARCH-001",
    ("src/core/tests", "REQ-002"): "REQ-ARCH-002",
    ("src/core/tests", "REQ-003"): "REQ-ARCH-003",
    ("src/core/tests", "REQ-004"): "REQ-ARCH-004",
    ("src/core/tests", "REQ-005"): "REQ-ARCH-005",
    ("src/core/tests", "REQ-006"): "REQ-ARCH-006",
    ("src/core/tests", "REQ-007"): "REQ-ARCH-007",
    # Physics test suite disambiguation (EP-4)
    ("tests/physics", "REQ-001"): "REQ-PHY-001",
    ("tests/physics", "REQ-002"): "REQ-PHY-002",
    ("tests/physics", "REQ-003"): "REQ-PHY-003",
    ("tests/physics", "REQ-004"): "REQ-PHY-004",
    ("tests/physics", "REQ-005"): "REQ-PHY-005",
    ("tests/physics", "REQ-006"): "REQ-PHY-006",
    ("tests/physics", "REQ-007"): "REQ-PHY-007",
    ("tests/physics", "REQ-008"): "REQ-PHY-008",
    ("tests/physics", "REQ-009"): "REQ-PHY-009",
    ("tests/physics", "REQ-010"): "REQ-PHY-010",
    ("tests/physics", "REQ-011"): "REQ-PHY-011",
    # Packaging legacy IDs (EP-2)
    (None, "REQ-PKG-01"): "REQ-PKG-001",
    (None, "REQ-PKG-02"): "REQ-PKG-002",
    (None, "REQ-PKG-03"): "REQ-PKG-003",
    (None, "REQ-PKG-04"): "REQ-PKG-004",
    (None, "REQ-PKG-05"): "REQ-PKG-005",
    (None, "REQ-PKG-06"): "REQ-PKG-006",
    (None, "REQ-PKG-07"): "REQ-PKG-007",
    (None, "REQ-PKG-08"): "REQ-PKG-008",
    (None, "REQ-PKG-09"): "REQ-PKG-009",
    # CI legacy IDs (EP-3)
    (None, "REQ-ARC-001"): "REQ-CI-001",
    (None, "REQ-ARC-002"): "REQ-CI-002",
    (None, "REQ-ARC-003"): "REQ-CI-003",
    (None, "REQ-ARC-004"): "REQ-CI-004",
    (None, "REQ-ARC-005"): "REQ-CI-005",
    (None, "REQ-ARC-006"): "REQ-CI-006",
    (None, "REQ-ARC-007"): "REQ-CI-007",
    (None, "REQ-ARC-008"): "REQ-CI-008",
    (None, "REQ-ARC-009"): "REQ-CI-009",
    (None, "REQ-ARC-010"): "REQ-CI-010",
}


def canonicalize_req(req_id: str, source_path: str = "") -> str:
    """Map legacy or un-namespaced requirement IDs to canonical namespaced IDs."""
    normalized_path = source_path.replace("\\", "/")
    # Check path-prefixed rules first
    for (prefix, leg_id), canonical in LEGACY_REQ_ALIASES.items():
        if leg_id == req_id and prefix is not None and normalized_path.startswith(prefix):
            return canonical
    # Check fallback rules with no prefix constraint
    for (prefix, leg_id), canonical in LEGACY_REQ_ALIASES.items():
        if leg_id == req_id and prefix is None:
            return canonical
    return req_id


@dataclasses.dataclass
class Requirement:
    """Represents a requirement parsed from an Enhancement Proposal."""

    req_id: str
    summary: str
    lifecycle_status: str
    invariant: str
    ep_file: str
    ep_num: int
    is_implemented: bool

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


@dataclasses.dataclass
class TestMapping:
    """Represents a link between a requirement and a test case."""

    req_id: str
    raw_id: str
    source_file: str
    line_number: int
    framework: str
    test_suite: str
    test_name: str
    full_symbol: str

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


def shield_math_pipes(text: str) -> str:
    """Shield Markdown table pipe delimiters occurring inside LaTeX math expressions."""
    def repl(m: re.Match[str]) -> str:
        return m.group(0).replace("|", "\x00PIPE\x00")

    text = re.sub(r"\$\$.*?\$\$", repl, text, flags=re.DOTALL)
    text = re.sub(r"\$.*?\$", repl, text)
    text = re.sub(r"\\\(.*?\\\)", repl, text)
    text = re.sub(r"\\\[.*?\\\]", repl, text, flags=re.DOTALL)
    return text


def restore_math_pipes(text: str) -> str:
    """Restore table pipe characters previously shielded inside math expressions."""
    return text.replace("\x00PIPE\x00", "|")


def parse_ep_file(ep_path: Path, repo_root: Path) -> tuple[list[Requirement], list[str]]:
    """Parse requirements from Section 3 of an Enhancement Proposal Markdown document."""
    requirements: list[Requirement] = []
    errors: list[str] = []

    rel_path = ep_path.relative_to(repo_root).as_posix()
    ep_match = re.search(r"ep-(\d+)", ep_path.name)
    ep_num = int(ep_match.group(1)) if ep_match else 0

    try:
        content = ep_path.read_text(encoding="utf-8")
    except Exception as exc:
        return [], [f"{rel_path}: Failed to read file: {exc}"]

    # Extract Section 3
    section3_match = re.search(r"(?m)^## 3\..*?$(.*?)(?=^## |\Z)", content, flags=re.DOTALL)
    if not section3_match:
        return [], [f"{rel_path}: Section 3 heading ('## 3.') not found"]

    section3_text = section3_match.group(1)
    lines = section3_text.splitlines()

    in_table = False
    headers: list[str] = []
    req_col = 0
    summary_col = 1
    lifecycle_col = 2
    invariant_col = 3

    for line_idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped.startswith("|") or not stripped.endswith("|"):
            in_table = False
            continue

        shielded = shield_math_pipes(stripped)
        raw_cells = [restore_math_pipes(c.strip()) for c in shielded.strip("|").split("|")]

        if not in_table:
            # Check for requirement table header
            norm_header = [c.lower() for c in raw_cells]
            has_req = any("requirement" in h or "req" in h for h in norm_header)
            has_summary = any("summary" in h or "name" in h or "title" in h for h in norm_header)
            if len(norm_header) >= 3 and has_req and has_summary:
                in_table = True
                headers = norm_header
                # Identify column positions
                for i, h in enumerate(headers):
                    if "requirement id" in h or "req id" in h or h == "requirement":
                        req_col = i
                    elif "summary" in h or "name" in h or "title" in h:
                        summary_col = i
                    elif "lifecycle" in h or "status" in h:
                        lifecycle_col = i
                    elif "invariant" in h or "definition" in h:
                        invariant_col = i
                continue
            continue

        # Skip separator line (e.g. | :--- | :--- |)
        if all(re.match(r"^:?-+:?$", c) for c in raw_cells):
            continue

        if len(raw_cells) <= max(req_col, summary_col):
            continue

        raw_req_cell = raw_cells[req_col]
        id_match = re.search(r"REQ-[A-Za-z0-9_-]+", raw_req_cell)
        if not id_match:
            continue

        raw_id = id_match.group(0)
        canonical_id = canonicalize_req(raw_id, rel_path)

        # Validate canonical ID against global naming standard REQ-TRC-001
        if not REQ_ID_PATTERN.match(canonical_id):
            errors.append(
                f"{rel_path}:{line_idx}: Invalid requirement ID format '{canonical_id}'. "
                f"Must match {REQ_ID_PATTERN.pattern}"
            )

        summary = raw_cells[summary_col].replace("**", "").strip()

        # Handle 3-column vs 4-column schemas
        if len(raw_cells) > lifecycle_col and lifecycle_col != invariant_col:
            lifecycle = raw_cells[lifecycle_col].strip()
        else:
            lifecycle = "Phase 1 (Planned)"

        if len(raw_cells) > invariant_col:
            invariant = raw_cells[invariant_col].strip()
        elif len(raw_cells) > 2:
            invariant = raw_cells[-1].strip()
        else:
            invariant = ""

        is_implemented = "(implemented)" in lifecycle.lower()

        requirements.append(
            Requirement(
                req_id=canonical_id,
                summary=summary,
                lifecycle_status=lifecycle,
                invariant=invariant,
                ep_file=rel_path,
                ep_num=ep_num,
                is_implemented=is_implemented,
            )
        )

    return requirements, errors


def extract_cpp_tests(file_path: Path, repo_root: Path) -> list[TestMapping]:
    """Extract requirement tags and test symbols from GoogleTest and Catch2 C++ files."""
    mappings: list[TestMapping] = []
    rel_path = file_path.relative_to(repo_root).as_posix()

    try:
        content = file_path.read_text(encoding="utf-8")
    except Exception:
        return []

    lines = content.splitlines()
    pending_tags: list[tuple[str, int]] = []
    non_comment_lines = 0

    gtest_pattern = re.compile(
        r"^\s*(TEST|TEST_F|TEST_P|TYPED_TEST|TYPED_TEST_P)\s*\(\s*([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*\)"
    )
    catch2_pattern = re.compile(
        r'^\s*TEST_CASE\s*\(\s*"([^"]+)"(?:\s*,\s*"([^"]*)")?\s*\)'
    )
    catch2_section_pattern = re.compile(
        r'^\s*SECTION\s*\(\s*"([^"]+)"(?:\s*,\s*"([^"]*)")?\s*\)'
    )

    for line_idx, line in enumerate(lines, start=1):
        stripped = line.strip()

        # Extract tags from line comments
        comment_match = re.search(r"//\s*(.*)", line)
        block_match = re.search(r"/\*\s*(.*?)\s*\*/", line)
        comment_text = ""
        if comment_match:
            comment_text = comment_match.group(1)
        elif block_match:
            comment_text = block_match.group(1)

        if comment_text:
            tags = REQ_TAG_SEARCH_PATTERN.findall(comment_text)
            for tag in tags:
                pending_tags.append((tag, line_idx))
            non_comment_lines = 0
            continue

        if not stripped:
            continue

        # Check for GoogleTest macro
        gm = gtest_pattern.match(stripped)
        if gm:
            suite_name = gm.group(2)
            test_name = gm.group(3)
            full_symbol = f"{suite_name}.{test_name}"

            for raw_tag, tag_line in pending_tags:
                can_id = canonicalize_req(raw_tag, rel_path)
                mappings.append(
                    TestMapping(
                        req_id=can_id,
                        raw_id=raw_tag,
                        source_file=rel_path,
                        line_number=tag_line,
                        framework="gtest",
                        test_suite=suite_name,
                        test_name=test_name,
                        full_symbol=full_symbol,
                    )
                )
            pending_tags = []
            non_comment_lines = 0
            continue

        # Check for Catch2 TEST_CASE
        cm = catch2_pattern.match(stripped)
        if cm:
            test_desc = cm.group(1)
            tag_str = cm.group(2) or ""
            embedded_tags = REQ_TAG_SEARCH_PATTERN.findall(tag_str)
            all_tags = pending_tags + [(t, line_idx) for t in embedded_tags]

            for raw_tag, tag_line in all_tags:
                can_id = canonicalize_req(raw_tag, rel_path)
                mappings.append(
                    TestMapping(
                        req_id=can_id,
                        raw_id=raw_tag,
                        source_file=rel_path,
                        line_number=tag_line,
                        framework="catch2",
                        test_suite=test_desc,
                        test_name=test_desc,
                        full_symbol=test_desc,
                    )
                )
            pending_tags = []
            non_comment_lines = 0
            continue

        # Check for Catch2 SECTION
        sm = catch2_section_pattern.match(stripped)
        if sm:
            section_desc = sm.group(1)
            tag_str = sm.group(2) or ""
            embedded_tags = REQ_TAG_SEARCH_PATTERN.findall(tag_str)
            all_tags = pending_tags + [(t, line_idx) for t in embedded_tags]

            for raw_tag, tag_line in all_tags:
                can_id = canonicalize_req(raw_tag, rel_path)
                mappings.append(
                    TestMapping(
                        req_id=can_id,
                        raw_id=raw_tag,
                        source_file=rel_path,
                        line_number=tag_line,
                        framework="catch2",
                        test_suite=section_desc,
                        test_name=section_desc,
                        full_symbol=section_desc,
                    )
                )
            pending_tags = []
            non_comment_lines = 0
            continue

        # Reset pending tags if non-comment lines accumulate without test macro
        non_comment_lines += 1
        if non_comment_lines > 3:
            pending_tags = []

    return mappings


def extract_python_tests(file_path: Path, repo_root: Path) -> list[TestMapping]:
    """Extract requirement tags and test symbols from Pytest files using Python AST."""
    mappings: list[TestMapping] = []
    rel_path = file_path.relative_to(repo_root).as_posix()

    try:
        content = file_path.read_text(encoding="utf-8")
        tree = ast.parse(content, filename=str(file_path))
    except Exception:
        return []

    lines = content.splitlines()

    # Module name representation
    mod_name = file_path.stem

    class TestVisitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.class_name: str | None = None

        def visit_ClassDef(self, node: ast.ClassDef) -> None:
            prev_class = self.class_name
            self.class_name = node.name
            self.generic_visit(node)
            self.class_name = prev_class

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            self._process_func(node)
            self.generic_visit(node)

        def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
            self._process_func(node)
            self.generic_visit(node)

        def _process_func(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
            found_tags: list[tuple[str, int]] = []

            # 1. Inspect Pytest markers: @pytest.mark.req("REQ-...")
            for decorator in node.decorator_list:
                call_node: ast.Call | None = None
                if isinstance(decorator, ast.Call):
                    call_node = decorator

                if call_node:
                    func_expr = call_node.func
                    is_req_marker = False
                    if isinstance(func_expr, ast.Attribute) and func_expr.attr == "req":
                        is_req_marker = True
                    elif isinstance(func_expr, ast.Name) and func_expr.id == "req":
                        is_req_marker = True

                    if is_req_marker:
                        for arg in call_node.args:
                            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                                found_tags.append((arg.value, node.lineno))

            # 2. Inspect docstring
            docstring = ast.get_docstring(node)
            if docstring:
                for match in REQ_TAG_SEARCH_PATTERN.finditer(docstring):
                    found_tags.append((match.group(1), node.lineno))

            # 3. Inspect preceding comments
            start_line = node.lineno - 1
            idx = start_line - 1
            while idx >= 0:
                line_str = lines[idx].strip()
                if line_str.startswith("#"):
                    for m in REQ_TAG_SEARCH_PATTERN.finditer(line_str):
                        found_tags.append((m.group(1), idx + 1))
                    idx -= 1
                else:
                    break

            if found_tags:
                test_name = node.name
                if self.class_name:
                    suite_name = f"{mod_name}.{self.class_name}"
                    full_symbol = f"{mod_name}.{self.class_name}.{test_name}"
                else:
                    suite_name = mod_name
                    full_symbol = f"{mod_name}.{test_name}"

                for raw_tag, tag_line in found_tags:
                    can_id = canonicalize_req(raw_tag, rel_path)
                    mappings.append(
                        TestMapping(
                            req_id=can_id,
                            raw_id=raw_tag,
                            source_file=rel_path,
                            line_number=tag_line,
                            framework="pytest",
                            test_suite=suite_name,
                            test_name=test_name,
                            full_symbol=full_symbol,
                        )
                    )

    visitor = TestVisitor()
    visitor.visit(tree)
    return mappings


def collect_all_requirements(
    eps_dir: Path, repo_root: Path
) -> tuple[dict[str, Requirement], list[str]]:
    """Scan all Enhancement Proposals and collect canonical requirements."""
    req_map: dict[str, Requirement] = {}
    errors: list[str] = []

    if not eps_dir.exists():
        return {}, [f"EP directory does not exist: {eps_dir}"]

    ep_files = sorted(eps_dir.glob("ep-*.md"))
    for ep_file in ep_files:
        if "template" in ep_file.name.lower():
            continue
        parsed_requirements, ep_errs = parse_ep_file(ep_file, repo_root)
        errors.extend(ep_errs)
        for req in parsed_requirements:
            if req.req_id in req_map:
                existing = req_map[req.req_id]
                errors.append(
                    f"Duplicate requirement ID '{req.req_id}' found in '{req.ep_file}' "
                    f"(previously defined in '{existing.ep_file}')"
                )
            else:
                req_map[req.req_id] = req

    return req_map, errors


def collect_all_tests(repo_root: Path) -> list[TestMapping]:
    """Scan the repository test suites and return all requirement mappings."""
    all_mappings: list[TestMapping] = []
    search_dirs = [repo_root / "src", repo_root / "tests", repo_root / "scripts" / "tests"]

    for search_dir in search_dirs:
        if not search_dir.exists():
            continue

        # C++ files
        for ext in ("*.cpp", "*.hpp", "*.h", "*.cc"):
            for cpp_file in search_dir.rglob(ext):
                all_mappings.extend(extract_cpp_tests(cpp_file, repo_root))

        # Python files
        for py_file in search_dir.rglob("*.py"):
            all_mappings.extend(extract_python_tests(py_file, repo_root))

    return all_mappings


def verify_manifest(manifest_path: Path, repo_root: Path) -> tuple[bool, list[str]]:
    """Verify sealed qualification acceptance test integrity using normalised LF line endings."""
    if not manifest_path.exists():
        return False, [f"Manifest file not found: {manifest_path}"]

    errors: list[str] = []
    lines = manifest_path.read_text(encoding="utf-8").splitlines()

    verified_count = 0
    for line_idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        parts = stripped.split(None, 1)
        if len(parts) != 2:
            errors.append(f"{manifest_path}:{line_idx}: Invalid manifest line format. Expected: '<hash> <filepath>'")
            continue

        expected_hash, rel_target = parts[0].strip(), parts[1].strip()
        target_path = repo_root / rel_target

        if not target_path.exists():
            errors.append(f"Sealed qualification file missing: {rel_target}")
            continue

        try:
            raw_bytes = target_path.read_bytes()
            # Normalise Windows CRLF to LF across OS runners
            normalised = raw_bytes.replace(b"\r\n", b"\n")
            actual_hash = hashlib.sha256(normalised).hexdigest()

            if actual_hash != expected_hash:
                errors.append(
                    f"Integrity violation in sealed file '{rel_target}': "
                    f"expected {expected_hash}, calculated {actual_hash}"
                )
            else:
                verified_count += 1
        except Exception as exc:
            errors.append(f"Failed to read sealed file '{rel_target}': {exc}")

    success = len(errors) == 0
    return success, errors


@dataclasses.dataclass
class JUnitResult:
    """Parsed JUnit XML test case result."""

    name: str
    classname: str
    status: str  # "PASSED", "FAILED", "ERROR", "SKIPPED"
    message: str = ""


def parse_junit_xml(xml_path: Path) -> dict[str, JUnitResult]:
    """Parse a JUnit XML file generated by CTest or Pytest."""
    results: dict[str, JUnitResult] = {}
    if not xml_path.exists():
        return results

    try:
        tree = ET.parse(xml_path)
    except Exception:
        return results

    root = tree.getroot()
    # Test cases can be directly under root or within <testsuite> elements
    for tc in root.iter("testcase"):
        name = tc.attrib.get("name", "")
        classname = tc.attrib.get("classname", "")

        status = "PASSED"
        msg = ""

        # Check child elements in order of severity
        err_node = tc.find("error")
        fail_node = tc.find("failure")
        skip_node = tc.find("skipped")

        if err_node is not None:
            status = "ERROR"
            msg = err_node.attrib.get("message", "") or (err_node.text or "").strip()
        elif fail_node is not None:
            status = "FAILED"
            msg = fail_node.attrib.get("message", "") or (fail_node.text or "").strip()
        elif skip_node is not None:
            status = "SKIPPED"
            msg = skip_node.attrib.get("message", "") or (skip_node.text or "").strip()

        # Key by multiple symbol forms to allow flexible matching
        key_full = f"{classname}.{name}" if classname else name
        res = JUnitResult(name=name, classname=classname, status=status, message=msg)
        results[key_full] = res
        results[name] = res

    return results


def match_test_result(
    mapping: TestMapping, xml_results: dict[str, JUnitResult]
) -> JUnitResult | None:
    """Match a static TestMapping to a dynamically executed JUnit test result."""
    # 1. Exact full_symbol match (e.g. DiracOperatorTests.NoErrorsWhenConstruction)
    if mapping.full_symbol in xml_results:
        return xml_results[mapping.full_symbol]

    # 2. Match by test_name (e.g. NoErrorsWhenConstruction)
    if mapping.test_name in xml_results:
        return xml_results[mapping.test_name]

    # 3. Match by suffix
    for key, res in xml_results.items():
        if key.endswith(f".{mapping.test_name}"):
            return res

    return None


def audit_repository(repo_root: Path, eps_dir: Path) -> int:
    """Audit the repository for duplicate IDs, orphaned tests, and specification coverage."""
    req_map, ep_errors = collect_all_requirements(eps_dir, repo_root)
    all_tests = collect_all_tests(repo_root)

    print("=" * 70)
    print(" RFL Specification Traceability & Governance Audit")
    print("=" * 70)

    # 1. Report EP parsing errors
    if ep_errors:
        print(f"\n❌ Found {len(ep_errors)} specification validation error(s):")
        for err in ep_errors:
            print(f"  • {err}")
    else:
        print(f"✅ Verified {len(req_map)} requirement IDs across Enhancement Proposals.")

    # 2. Check for orphaned tests
    orphaned = [t for t in all_tests if t.req_id not in req_map]
    if orphaned:
        print(f"\n⚠️  Found {len(orphaned)} orphaned test mapping(s) referencing unknown requirements:")
        for t in orphaned:
            print(f"  • [{t.req_id}] in {t.source_file}:{t.line_number} ({t.full_symbol})")

    # 3. Compute requirement coverage
    implemented_requirements = [r for r in req_map.values() if r.is_implemented]
    planned_requirements = [r for r in req_map.values() if not r.is_implemented]

    req_to_tests: dict[str, list[TestMapping]] = {r_id: [] for r_id in req_map}
    for t in all_tests:
        if t.req_id in req_to_tests:
            req_to_tests[t.req_id].append(t)

    covered_implemented = [r for r in implemented_requirements if req_to_tests[r.req_id]]
    missing_implemented = [r for r in implemented_requirements if not req_to_tests[r.req_id]]

    print("\n📊 Traceability Coverage Statistics:")
    print(f"  • Total Requirements:   {len(req_map)}")
    print(f"  • Implemented:          {len(implemented_requirements)}")
    print(f"  • Planned / In-Flight:  {len(planned_requirements)}")
    print(f"  • Tagged Test Mappings: {len(all_tests)}")
    if implemented_requirements:
        cov_pct = (len(covered_implemented) / len(implemented_requirements)) * 100.0
        print(f"  • Implemented Coverage: {len(covered_implemented)}/{len(implemented_requirements)} ({cov_pct:.1f}%)")

    if missing_implemented:
        print(f"\n⚠️  Implemented requirements lacking mapped test cases ({len(missing_implemented)}):")
        for r in missing_implemented:
            print(f"  • [{r.req_id}] {r.summary} ({r.ep_file})")

    print("\n" + "=" * 70)

    # Exit code: 0 if valid IDs and zero EP errors; 1 if errors found
    if ep_errors:
        return 1
    return 0


def query_requirement(req_id: str, repo_root: Path, eps_dir: Path) -> int:
    """Display detailed requirement metadata and mapped test cases."""
    req_map, _ = collect_all_requirements(eps_dir, repo_root)
    all_tests = collect_all_tests(repo_root)

    canonical_id = canonicalize_req(req_id)
    if canonical_id not in req_map:
        print(f"❌ Requirement '{req_id}' (canonical: '{canonical_id}') not found in any EP.")
        return 1

    req = req_map[canonical_id]
    mapped_tests = [t for t in all_tests if t.req_id == canonical_id]

    print("=" * 70)
    print(f" Requirement: {req.req_id}")
    print("=" * 70)
    print(f"• Summary:          {req.summary}")
    print(f"• Lifecycle Status: {req.lifecycle_status}")
    print(f"• Specification:    {req.ep_file} (EP-{req.ep_num})")
    print(f"• Invariant:        {req.invariant}")
    print(f"\n• Mapped Tests ({len(mapped_tests)}):")
    if mapped_tests:
        for t in mapped_tests:
            print(f"  - [{t.framework}] {t.full_symbol} ({t.source_file}:{t.line_number})")
    else:
        print("  - None")
    print("=" * 70)
    return 0


def generate_traceability_matrix(
    output_path: Path, repo_root: Path, eps_dir: Path, xml_paths: list[Path] | None = None
) -> None:
    """Generate a GitHub Flavoured Markdown traceability matrix document."""
    req_map, _ = collect_all_requirements(eps_dir, repo_root)
    all_tests = collect_all_tests(repo_root)

    xml_results: dict[str, JUnitResult] = {}
    if xml_paths:
        for p in xml_paths:
            xml_results.update(parse_junit_xml(p))

    req_to_tests: dict[str, list[TestMapping]] = {r_id: [] for r_id in req_map}
    for t in all_tests:
        if t.req_id in req_to_tests:
            req_to_tests[t.req_id].append(t)

    lines: list[str] = [
        "# RFL Requirements Traceability Matrix",
        "",
        "> Automated living compliance matrix generated by `scripts/verify_traceability.py`.",
        "",
        "| Requirement ID | Summary | EP | Lifecycle Status | Mapped Tests | Dynamic Execution |",
        "| :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for req_id in sorted(req_map.keys()):
        req = req_map[req_id]
        tests = req_to_tests[req_id]
        test_count = len(tests)

        if tests:
            test_links = "<br>".join(f"`{t.full_symbol}`" for t in tests[:3])
            if len(tests) > 3:
                test_links += f"<br>*(+{len(tests) - 3} more)*"
        else:
            test_links = "—"

        dynamic_status = "⚪ Not Tested"
        if xml_results and tests:
            statuses = [match_test_result(t, xml_results) for t in tests]
            valid_statuses = [s for s in statuses if s is not None]
            if not valid_statuses:
                dynamic_status = "⚪ Missing XML"
            elif any(s.status == "ERROR" for s in valid_statuses):
                dynamic_status = "💥 Error"
            elif any(s.status == "FAILED" for s in valid_statuses):
                dynamic_status = "❌ Failed"
            elif any(s.status == "SKIPPED" for s in valid_statuses):
                dynamic_status = "⚠️ Skipped"
            else:
                dynamic_status = "✅ Passing"
        elif req.is_implemented:
            dynamic_status = "⏳ Implemented"
        else:
            dynamic_status = "💡 Planned"

        ep_link = f"EP-{req.ep_num}"
        target_abs = repo_root / req.ep_file
        rel_target = os.path.relpath(target_abs, output_path.parent).replace("\\", "/")
        lines.append(
            f"| **{req.req_id}** | {req.summary} | [{ep_link}]({rel_target}) | "
            f"{req.lifecycle_status} | {test_links} | {dynamic_status} |"
        )

    lines.append("")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"✅ Generated traceability matrix at: {output_path}")


def dump_json(
    repo_root: Path, eps_dir: Path, xml_paths: list[Path] | None = None
) -> None:
    """Dump machine-readable JSON representation for AI agents and tooling."""
    req_map, _ = collect_all_requirements(eps_dir, repo_root)
    all_tests = collect_all_tests(repo_root)

    xml_results: dict[str, JUnitResult] = {}
    if xml_paths:
        for p in xml_paths:
            xml_results.update(parse_junit_xml(p))

    req_data = []
    for r_id, req in sorted(req_map.items()):
        tests = [t.to_dict() for t in all_tests if t.req_id == r_id]
        req_data.append({
            **req.to_dict(),
            "tests": tests,
            "test_count": len(tests),
        })

    payload = {
        "requirements": req_data,
        "total_requirements": len(req_map),
        "total_test_mappings": len(all_tests),
    }
    print(json.dumps(payload, indent=2))


def check_pull_request(
    base_ref: str,
    repo_root: Path,
    eps_dir: Path,
    xml_paths: list[Path] | None = None,
    red_anchor: bool = False,
) -> int:
    """Verify Git trailers and test execution for a feature pull request."""
    req_map, _ = collect_all_requirements(eps_dir, repo_root)
    all_tests = collect_all_tests(repo_root)

    # Extract Git trailers: "Implements: REQ-..."
    cmd = ["git", "log", f"{base_ref}..HEAD", "--pretty=format:%B---COMMIT---"]
    try:
        process_result = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True, check=True)
        commit_texts = process_result.stdout.split("---COMMIT---")
    except Exception:
        # Fallback to inspecting HEAD commit alone
        process_result = subprocess.run(["git", "log", "-1", "--pretty=format:%B"], cwd=repo_root, capture_output=True, text=True)
        commit_texts = [process_result.stdout]

    claimed_raw_ids: set[str] = set()
    trailer_regex = re.compile(r"(?im)^Implements:\s*(.+)$")
    for txt in commit_texts:
        for m in trailer_regex.finditer(txt):
            trailer_val = m.group(1).strip()
            for token in REQ_TAG_SEARCH_PATTERN.findall(trailer_val):
                claimed_raw_ids.add(token)

    claimed_requirements = [canonicalize_req(r) for r in claimed_raw_ids]

    print("=" * 70)
    print(" RFL Pull Request Traceability Qualification Gate")
    print("=" * 70)

    if not claimed_requirements:
        print("ℹ️  No 'Implements: REQ-...' trailer found.")
        print("   Treated as internal refactoring (non-blocking informational pass).")
        print("=" * 70)
        return 0

    print(f"🔍 Claimed Requirement Delivery: {', '.join(claimed_requirements)}")

    errors: list[str] = []
    xml_results: dict[str, JUnitResult] = {}
    if xml_paths:
        for p in xml_paths:
            xml_results.update(parse_junit_xml(p))

    summary_rows: list[str] = [
        "| Requirement ID | Summary | Mapped Tests | Red Check | Green Check | Status |",
        "| :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for req_id in claimed_requirements:
        if req_id not in req_map:
            errors.append(f"Claimed requirement '{req_id}' does not exist in any EP.")
            summary_rows.append(f"| **{req_id}** | *Unknown* | — | — | — | ❌ Missing Spec |")
            continue

        req = req_map[req_id]
        tests = [t for t in all_tests if t.req_id == req_id]

        if not tests:
            errors.append(f"Requirement '{req_id}' has no mapped unit or acceptance test cases.")
            summary_rows.append(f"| **{req_id}** | {req.summary} | 0 tests | — | — | ❌ No Tests |")
            continue

        test_count_str = f"{len(tests)} test(s)"

        if red_anchor and xml_results:
            # Red Anchor verification: claimed tests must fail with <failure>, never <error>
            test_results = [match_test_result(t, xml_results) for t in tests]
            has_error = any(r is not None and r.status == "ERROR" for r in test_results)
            has_failure = any(r is not None and r.status == "FAILED" for r in test_results)
            if has_error:
                errors.append(f"Red Anchor check failed for '{req_id}': produced uncaught <error>.")
                summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | 💥 Error | — | ❌ Invalid Red |")
            elif not has_failure:
                errors.append(f"Red Anchor check failed for '{req_id}': expected <failure> on stubs, but test passed.")
                summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | ⚪ Green | — | ❌ Missing Red |")
            else:
                summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | 🔴 Failure | — | ✅ Valid Red |")
        elif xml_results:
            # Green Anchor verification: claimed tests must pass with zero failures
            test_results = [match_test_result(t, xml_results) for t in tests]
            missing_exec = any(r is None for r in test_results)
            has_error = any(r is not None and r.status == "ERROR" for r in test_results)
            has_failure = any(r is not None and r.status == "FAILED" for r in test_results)
            has_skipped = any(r is not None and r.status == "SKIPPED" for r in test_results)

            if missing_exec:
                errors.append(f"Target test for '{req_id}' did not execute in JUnit test run.")
                summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | — | ⚪ Missing | ❌ Unverified |")
            elif has_error or has_failure or has_skipped:
                errors.append(f"Target test execution failed for '{req_id}'.")
                summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | — | ❌ Failed | ❌ Green Failure |")
            else:
                summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | — | 🟢 Passed | ✅ Verified |")
        else:
            summary_rows.append(f"| **{req_id}** | {req.summary} | {test_count_str} | — | — | ⏳ Static Only |")

    # Write Step Summary if running in GitHub Actions
    step_summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary_path:
        with open(step_summary_path, "a", encoding="utf-8") as f:
            f.write("\n### 📋 Requirements Traceability Verification Summary\n\n")
            f.write("\n".join(summary_rows))
            f.write("\n\n")

    if errors:
        print("\n❌ Pull request verification failed:")
        for err in errors:
            print(f"  • {err}")
        print("=" * 70)
        return 1

    print("\n✅ All claimed requirements statically and dynamically verified.")
    print("=" * 70)
    return 0


def assert_zero_copy_fortran_aligned(array: Any, expected_ptr: int | None = None) -> bool:
    """Verify that a NumPy array wrapping an Armadillo matrix is Fortran-contiguous.

    Prevents silent matrix transposition and complex conjugation (H -> H*)
    of Hermitian operators across the C++ and Python interface (REQ-TRC-008).
    """
    if not hasattr(array, "flags"):
        raise AssertionError("Object is not a NumPy array")
    if not array.flags["F_CONTIGUOUS"]:
        raise AssertionError("Matrix is not Fortran-contiguous (column-major order required for Armadillo)")
    if expected_ptr is not None:
        actual_ptr = array.ctypes.data
        if actual_ptr != expected_ptr:
            raise AssertionError(f"Zero-copy violation: memory pointer {actual_ptr} != {expected_ptr}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RFL Requirements Traceability & Governance Verification Engine (EP-5)"
    )
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parent.parent,
                        help="Repository root directory")
    parser.add_argument("--eps-dir", type=Path, default=None,
                        help="Directory containing Enhancement Proposals (default: <repo-root>/docs/eps)")
    parser.add_argument("--audit", action="store_true",
                        help="Audit repository for duplicate IDs, orphaned tests, and coverage")
    parser.add_argument("--query", type=str, metavar="REQ_ID",
                        help="Query requirement summary, status, invariant, and mapped test cases")
    parser.add_argument("--check-pr", action="store_true",
                        help="Verify current PR branch Git trailers and mapped test execution")
    parser.add_argument("--base", type=str, default="origin/main",
                        help="Git base reference for PR verification (default: origin/main)")
    parser.add_argument("--verify-xml", nargs="+", type=Path, metavar="XML_PATH",
                        help="Verify test execution against JUnit XML outputs")
    parser.add_argument("--verify-manifest", type=Path, metavar="MANIFEST_PATH",
                        help="Verify cryptographic SHA-256 integrity of sealed qualification tests")
    parser.add_argument("--generate-matrix", type=Path, metavar="OUTPUT_MD",
                        help="Generate Markdown traceability matrix document")
    parser.add_argument("--dump-json", action="store_true",
                        help="Dump machine-readable JSON representation for AI agents")
    parser.add_argument("--red-anchor", action="store_true",
                        help="Assert that claimed tests produce runtime <failure> on degenerate stubs")

    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    eps_dir = (args.eps_dir or (repo_root / "docs" / "eps")).resolve()

    if args.verify_manifest:
        success, errors = verify_manifest(args.verify_manifest, repo_root)
        if not success:
            print(f"❌ Manifest verification failed for {args.verify_manifest}:")
            for err in errors:
                print(f"  • {err}")
            return 1
        print(f"✅ Sealed qualification manifest verified successfully: {args.verify_manifest}")
        return 0

    if args.query:
        return query_requirement(args.query, repo_root, eps_dir)

    if args.generate_matrix:
        generate_traceability_matrix(args.generate_matrix, repo_root, eps_dir, args.verify_xml)
        return 0

    if args.dump_json:
        dump_json(repo_root, eps_dir, args.verify_xml)
        return 0

    if args.check_pr:
        return check_pull_request(args.base, repo_root, eps_dir, args.verify_xml, red_anchor=args.red_anchor)

    # Default action: run audit
    return audit_repository(repo_root, eps_dir)


if __name__ == "__main__":
    sys.exit(main())
