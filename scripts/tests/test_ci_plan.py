#!/usr/bin/env python3
"""
Unit tests for scripts/ci_plan.py (EP-5: Specification-Driven Agentic Workflow).

Verifies:
- C++ #include dependency graph and transitive header/source traversal
- pybind11 AST bridge mapping C++ headers to Python binding tests
- Two-Hop requirement impact inference from source changes and EP diffs
- Full-run escalation on root build config changes (CMakeLists.txt)
- Python dependency config impact (pyproject.toml -> all Python tests)
- Self-mapping of scripts/ci_plan.py to scripts/tests/test_ci_plan.py
- GitHub Actions outputs (GITHUB_OUTPUT), milestone verification (check-milestone),
  and gatekeeper evaluation (gate)
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPT_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import ci_plan  # noqa: E402


def req(*req_ids: str):
    """Decorator compatible with extract_python_tests (@req('REQ-...'))."""

    def decorator(obj):
        return obj

    return decorator


@req("REQ-TRC-005")
class TestCIPlannerAndDependencyGraph(unittest.TestCase):
    """Test C++ #include graph, pybind11 AST bridge, and Two-Hop requirement impact planner [REQ-TRC-005]."""

    def test_pybind11_symbol_map_extraction(self) -> None:
        symbol_map = ci_plan.parse_pybind11_symbol_map(REPO_ROOT)
        self.assertEqual(symbol_map.get("DiracOperator"), "src/core/DiracOperator.hpp")
        self.assertEqual(symbol_map.get("Action"), "src/core/BarrettGlaser/Action.hpp")
        self.assertEqual(symbol_map.get("Metropolis"), "src/core/BarrettGlaser/Metropolis.hpp")
        self.assertEqual(symbol_map.get("StdRng"), "src/core/StdRng.hpp")
        self.assertEqual(symbol_map.get("GslRng"), "src/core/StdRng.hpp")
        self.assertEqual(symbol_map.get("set_max_clifford_mode"), "src/core/Clifford.hpp")

    def test_cpp_component_isolation_hamiltonian(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["src/core/BarrettGlaser/Hamiltonian.cpp"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertTrue(plan.has_changes)
        self.assertFalse(plan.is_full_run)
        self.assertEqual(plan.affected_cpp_suites, ["HamiltonianTests"])
        self.assertEqual(plan.ctest_filter, r"^(HamiltonianTests)\.")
        # Hamiltonian is not wrapped in bindings.cpp, so Python tests are unaffected
        self.assertEqual(plan.affected_python_test_files, [])

    def test_two_hop_requirement_impact_from_source_change(self) -> None:
        # Touching Clifford.cpp (which contains zero REQ- tags) must infer REQ-ARCH-001,
        # REQ-ARCH-003, and REQ-ARCH-006 via dependent C++ tests and Python integration tests.
        plan = ci_plan.compute_plan_from_files(
            changed_files=["src/core/Clifford.cpp"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertIn("CliffordTests", plan.affected_cpp_suites)
        self.assertIn("DiracOperatorTests", plan.affected_cpp_suites)
        self.assertIn("REQ-ARCH-001", plan.impacted_requirements)
        self.assertIn("REQ-ARCH-003", plan.impacted_requirements)
        self.assertIn("REQ-ARCH-006", plan.impacted_requirements)
        self.assertIn("component-test", plan.requirement_triggers["REQ-ARCH-001"])
        self.assertIn("component-test", plan.requirement_triggers["REQ-ARCH-006"])
        # Python tests using rfl.DiracOperator / rfl.set_max_clifford_mode are also affected via pybind11 bridge
        self.assertIn(
            "src/python_bindings/tests/integration/test_boundary_safety.py",
            plan.affected_python_test_files,
        )

    def test_ep_spec_diff_triggers_mapped_tests_without_code_change(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["docs/eps/ep-1-core-architecture-modernisation.md"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
            ep_diff_requirements={"REQ-ARCH-001"},
        )
        self.assertIn("REQ-ARCH-001", plan.impacted_requirements)
        self.assertIn("ep-spec-diff", plan.requirement_triggers["REQ-ARCH-001"])
        self.assertIn("CliffordTests", plan.affected_cpp_suites)
        self.assertIn("DiracOperatorTests", plan.affected_cpp_suites)
        self.assertTrue(plan.builds.test_ubuntu)
        self.assertFalse(plan.builds.test_macos)
        self.assertFalse(plan.builds.test_windows)

    def test_docs_prose_only_change_skips_cpp_builds(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["README.md"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertTrue(plan.builds.validate_docs)
        self.assertFalse(plan.builds.lint_cpp)
        self.assertFalse(plan.builds.test_ubuntu)
        self.assertFalse(plan.builds.test_macos)
        self.assertFalse(plan.builds.test_windows)
        self.assertEqual(plan.required_checks, ["validate_docs"])

    def test_full_run_plan(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=[],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
            full_run=True,
        )
        self.assertTrue(plan.is_full_run)
        self.assertEqual(plan.ctest_filter, "")
        self.assertEqual(
            plan.required_checks,
            ["validate_docs", "lint", "test_ubuntu", "test_macos", "test_windows"],
        )

    def test_cmake_change_triggers_full_cross_platform_run(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["CMakeLists.txt"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertEqual(plan.ctest_filter, "")
        self.assertIn("CliffordTests", plan.affected_cpp_suites)
        self.assertIn("HamiltonianTests", plan.affected_cpp_suites)
        self.assertTrue(plan.builds.test_ubuntu)
        self.assertTrue(plan.builds.test_macos)
        self.assertTrue(plan.builds.test_windows)

    def test_bindings_cpp_change_triggers_python_binding_tests(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["src/python_bindings/bindings.cpp"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertTrue(plan.needs_python_bindings_build)
        self.assertIn(
            "src/python_bindings/tests/unit/test_api_surface.py",
            plan.affected_python_test_files,
        )
        self.assertIn("REQ-ARCH-006", plan.impacted_requirements)
        self.assertTrue(plan.builds.lint_cpp)
        self.assertTrue(plan.builds.test_ubuntu)

    def test_pyproject_toml_change_triggers_all_python_tests(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["pyproject.toml"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertIn("ApiContractTest", plan.affected_cpp_suites)
        self.assertTrue(plan.needs_python_bindings_build)
        self.assertIn("scripts/tests/test_ci_plan.py", plan.affected_python_test_files)
        self.assertIn("scripts/tests/test_verify_traceability.py", plan.affected_python_test_files)
        self.assertIn(
            "src/python_bindings/tests/unit/test_api_surface.py",
            plan.affected_python_test_files,
        )
        self.assertIn(
            "src/python_bindings/tests/integration/test_boundary_safety.py",
            plan.affected_python_test_files,
        )
        self.assertIn("REQ-TRC-005", plan.impacted_requirements)
        self.assertIn("REQ-ARCH-006", plan.impacted_requirements)
        self.assertTrue(plan.builds.validate_docs)
        self.assertFalse(plan.builds.lint_cpp)
        self.assertTrue(plan.builds.test_ubuntu)

    def test_ci_plan_script_change_triggers_test_ci_plan(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["scripts/ci_plan.py"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertIn("scripts/tests/test_ci_plan.py", plan.affected_python_test_files)
        self.assertIn("REQ-TRC-005", plan.impacted_requirements)
        self.assertTrue(plan.builds.validate_docs)
        self.assertFalse(plan.builds.test_ubuntu)

    def test_write_github_output(self) -> None:
        plan = ci_plan.compute_plan_from_files(
            changed_files=["src/core/BarrettGlaser/Hamiltonian.cpp"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        with tempfile.NamedTemporaryFile(mode="r+", encoding="utf-8") as tmp:
            old_env = os.environ.get("GITHUB_OUTPUT")
            os.environ["GITHUB_OUTPUT"] = tmp.name
            try:
                ci_plan.write_github_output(plan)
            finally:
                if old_env is None:
                    os.environ.pop("GITHUB_OUTPUT", None)
                else:
                    os.environ["GITHUB_OUTPUT"] = old_env
            tmp.seek(0)
            content = tmp.read()
        self.assertIn("has_changes=true", content)
        self.assertIn("test_ubuntu=true", content)
        self.assertIn("ctest_filter=^(HamiltonianTests)\\.", content)

    def test_cmd_gate_evaluation(self) -> None:
        old_needs = os.environ.get("NEEDS_JSON")
        try:
            # 1. All required checks pass
            os.environ["NEEDS_JSON"] = json.dumps(
                {
                    "plan": {
                        "result": "success",
                        "outputs": {"required_checks": '["validate_docs", "test_ubuntu"]'},
                    },
                    "check_milestone": {"result": "success"},
                    "validate_docs": {"result": "success"},
                    "test_ubuntu": {"result": "success"},
                    "test_macos": {"result": "skipped"},
                }
            )
            self.assertEqual(ci_plan.cmd_gate(), 0)

            # 2. Required check fails
            os.environ["NEEDS_JSON"] = json.dumps(
                {
                    "plan": {
                        "result": "success",
                        "outputs": {"required_checks": '["validate_docs", "test_ubuntu"]'},
                    },
                    "check_milestone": {"result": "success"},
                    "validate_docs": {"result": "success"},
                    "test_ubuntu": {"result": "failure"},
                }
            )
            self.assertEqual(ci_plan.cmd_gate(), 1)

            # 3. Milestone check fails
            os.environ["NEEDS_JSON"] = json.dumps(
                {
                    "plan": {
                        "result": "success",
                        "outputs": {"required_checks": '["validate_docs"]'},
                    },
                    "check_milestone": {"result": "failure"},
                    "validate_docs": {"result": "success"},
                }
            )
            self.assertEqual(ci_plan.cmd_gate(), 1)
        finally:
            if old_needs is None:
                os.environ.pop("NEEDS_JSON", None)
            else:
                os.environ["NEEDS_JSON"] = old_needs

    def test_cmd_check_milestone(self) -> None:
        old_event_name = os.environ.get("GITHUB_EVENT_NAME")
        old_event_path = os.environ.get("GITHUB_EVENT_PATH")
        old_repo = os.environ.get("GITHUB_REPOSITORY")
        try:
            # 1. Non-PR event skips milestone check
            os.environ["GITHUB_EVENT_NAME"] = "workflow_dispatch"
            self.assertEqual(ci_plan.cmd_check_milestone(), 0)

            # 2. PR event with milestone assigned in event payload passes
            os.environ["GITHUB_EVENT_NAME"] = "pull_request"
            os.environ["GITHUB_REPOSITORY"] = ""
            with tempfile.NamedTemporaryFile(mode="w", suffix=".json", encoding="utf-8") as tmp:
                json.dump(
                    {"pull_request": {"number": 84, "milestone": {"title": "v0.3.0"}}},
                    tmp,
                )
                tmp.flush()
                os.environ["GITHUB_EVENT_PATH"] = tmp.name
                self.assertEqual(ci_plan.cmd_check_milestone(), 0)

            # 3. PR event without milestone fails
            with tempfile.NamedTemporaryFile(mode="w", suffix=".json", encoding="utf-8") as tmp:
                json.dump(
                    {"pull_request": {"number": 84, "milestone": None}},
                    tmp,
                )
                tmp.flush()
                os.environ["GITHUB_EVENT_PATH"] = tmp.name
                self.assertEqual(ci_plan.cmd_check_milestone(), 1)
        finally:
            for k, v in (
                ("GITHUB_EVENT_NAME", old_event_name),
                ("GITHUB_EVENT_PATH", old_event_path),
                ("GITHUB_REPOSITORY", old_repo),
            ):
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v


if __name__ == "__main__":
    unittest.main()
