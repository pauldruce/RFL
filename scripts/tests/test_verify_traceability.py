#!/usr/bin/env python3
"""scripts/tests/test_verify_traceability.py

Unit test suite for the RFL Requirements Traceability and Governance Engine.
Verifies requirement ID validation, EP table parsing, domain-aware coverage classification,
C++ (GTest/Catch2) extraction, Python AST extraction (ignoring non-test helper functions),
JUnit XML evaluation, and shift-left hook integrity.

Conforms to EP-5 (REQ-TRC-001 through REQ-TRC-006).
Uses standard Python library only (runnable via unittest or pytest).
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import verify_traceability as vt


class TestRequirementIDValidation(unittest.TestCase):
    """Test requirement ID regex and namespacing rules [REQ-TRC-001]."""

    def test_valid_canonical_ids(self) -> None:
        valid_ids = [
            "REQ-ARCH-001",
            "REQ-PKG-009",
            "REQ-CI-010",
            "REQ-PHY-001",
            "REQ-TRC-006",
            "REQ-AI-042",
        ]
        for req_id in valid_ids:
            self.assertIsNotNone(
                vt.REQ_ID_PATTERN.match(req_id),
                f"Expected valid ID '{req_id}' to match pattern",
            )

    def test_invalid_canonical_ids(self) -> None:
        invalid_ids = [
            "REQ-001",          # Un-namespaced legacy format
            "REQ-A-001",        # Subsystem code too short (< 2 letters)
            "REQ-INVALID-001",  # Subsystem code too long (> 4 letters)
            "REQ-arch-001",     # Lowercase subsystem code
            "REQ-ARCH-1",       # Non-zero-padded number
            "REQ-ARCH-01",      # 2-digit number
            "ARCH-001",         # Missing REQ prefix
            "REQ-008b",         # Alphanumeric suffix
        ]
        for req_id in invalid_ids:
            self.assertIsNone(
                vt.REQ_ID_PATTERN.match(req_id),
                f"Expected invalid ID '{req_id}' to fail pattern match",
            )


class TestLegacyRequirementAliasing(unittest.TestCase):
    """Test legacy alias resolution and collision disambiguation [REQ-TRC-003]."""

    def test_path_prefixed_disambiguation(self) -> None:
        arch_id = vt.canonicalize_req("REQ-001", "src/core/tests/tDiracOperator.cpp")
        self.assertEqual(arch_id, "REQ-ARCH-001")

        delta_phy_id = vt.canonicalize_req("REQ-001", "src/core/tests/tDelta.cpp")
        self.assertEqual(delta_phy_id, "REQ-PHY-001")

        phy_id = vt.canonicalize_req("REQ-001", "tests/physics/test_invariants.py")
        self.assertEqual(phy_id, "REQ-PHY-001")

    def test_unambiguous_legacy_aliases(self) -> None:
        pkg_id = vt.canonicalize_req("REQ-PKG-01", "any/arbitrary/path.py")
        self.assertEqual(pkg_id, "REQ-PKG-001")

        ci_id = vt.canonicalize_req("REQ-ARC-001", "")
        self.assertEqual(ci_id, "REQ-CI-001")

        phy_011 = vt.canonicalize_req("REQ-011", "src/core/tests/tDelta.cpp")
        self.assertEqual(phy_011, "REQ-PHY-011")

    def test_already_canonical_ids_unchanged(self) -> None:
        modern_id = vt.canonicalize_req("REQ-ARCH-001", "src/core/tests/tDiracOperator.cpp")
        self.assertEqual(modern_id, "REQ-ARCH-001")


class TestMathPipeShielding(unittest.TestCase):
    """Test LaTeX math pipe character shielding to preserve Markdown table integrity [REQ-TRC-002]."""

    def test_single_dollar_math_shielding(self) -> None:
        text = "Formula $|S|$ and $\\lvert A \\rvert \\le B$ in text"
        shielded = vt.shield_math_pipes(text)
        self.assertNotIn("|", shielded)
        restored = vt.restore_math_pipes(shielded)
        self.assertEqual(restored, text)

    def test_double_dollar_math_shielding(self) -> None:
        text = "$$\\lVert D - D^\\dagger \\rVert_F \\le c_2 (p+q) \\epsilon_{\\mathrm{mach}} \\lVert D \\rVert_F$$"
        shielded = vt.shield_math_pipes(text)
        self.assertNotIn("|", shielded)
        restored = vt.restore_math_pipes(shielded)
        self.assertEqual(restored, text)


class TestEPTableParserAndDomainAwareness(unittest.TestCase):
    """Test Markdown requirement table extraction and domain classification [REQ-TRC-002]."""

    def test_parse_four_column_ep_table(self) -> None:
        ep_content = """# EP-99: Test Proposal

## 3. Research Workflows & Scientific Requirements

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Lifecycle Status | Operational & Physical Invariant |
| :--- | :--- | :--- | :--- |
| **REQ-ARCH-001** | **Hermiticity Preservation** | Phase 1 (Implemented) | Matrix variations $\\delta M_k$ must preserve $H_i^\\dagger = H_i$. |
| **REQ-ARCH-002** | **Analytic Variation** | Phase 2 (Planned) | Fast $\\mathcal{O}(N^3)$ local variation equals global recomputation. |
| **REQ-CI-001** | **Fast Dual-Platform PR Gate** | Phase 1 (Implemented) | PR gate finishes in under 120 seconds. |
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            ep_file = tmp_root / "ep-99-test.md"
            ep_file.write_text(ep_content, encoding="utf-8")

            parsed_requirements, errors = vt.parse_ep_file(ep_file, tmp_root)
            self.assertEqual(len(errors), 0, f"Errors: {errors}")
            self.assertEqual(len(parsed_requirements), 3)

            self.assertEqual(parsed_requirements[0].req_id, "REQ-ARCH-001")
            self.assertTrue(parsed_requirements[0].is_implemented)
            self.assertTrue(parsed_requirements[0].is_unit_testable)

            self.assertEqual(parsed_requirements[1].req_id, "REQ-ARCH-002")
            self.assertFalse(parsed_requirements[1].is_implemented)
            self.assertTrue(parsed_requirements[1].is_unit_testable)

            self.assertEqual(parsed_requirements[2].req_id, "REQ-CI-001")
            self.assertTrue(parsed_requirements[2].is_implemented)
            self.assertFalse(parsed_requirements[2].is_unit_testable)


class TestCppTestExtraction(unittest.TestCase):
    """Test requirement tag extraction from GoogleTest and Catch2 source files [REQ-TRC-004]."""

    def test_extract_googletest_comments(self) -> None:
        cpp_content = """
#include <gtest/gtest.h>

// [REQ-ARCH-001]
TEST(DiracOperatorTests, PreservesHermiticity) {
    EXPECT_TRUE(true);
}

// [REQ-ARCH-002] Multi-tag test
// [REQ-ARCH-003]
TEST_F(DiracFixture, EvaluatesVariation) {
    EXPECT_TRUE(true);
}
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            src_file = tmp_root / "tDirac.cpp"
            src_file.write_text(cpp_content, encoding="utf-8")

            mappings = vt.extract_cpp_tests(src_file, tmp_root)
            self.assertEqual(len(mappings), 3)

            req_symbols = {(m.req_id, m.full_symbol) for m in mappings}
            self.assertIn(("REQ-ARCH-001", "DiracOperatorTests.PreservesHermiticity"), req_symbols)
            self.assertIn(("REQ-ARCH-002", "DiracFixture.EvaluatesVariation"), req_symbols)
            self.assertIn(("REQ-ARCH-003", "DiracFixture.EvaluatesVariation"), req_symbols)

    def test_extract_catch2_bracketed_tags(self) -> None:
        cpp_content = """
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Dirac Operator Hermiticity", "[core][REQ-ARCH-001]") {
    SECTION("Sub-variation test", "[REQ-ARCH-002]") {
        REQUIRE(true);
    }
}
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            src_file = tmp_root / "tCatch.cpp"
            src_file.write_text(cpp_content, encoding="utf-8")

            mappings = vt.extract_cpp_tests(src_file, tmp_root)
            self.assertEqual(len(mappings), 2)
            req_ids = {m.req_id for m in mappings}
            self.assertIn("REQ-ARCH-001", req_ids)
            self.assertIn("REQ-ARCH-002", req_ids)


class TestPythonTestExtraction(unittest.TestCase):
    """Test requirement tag extraction from Pytest files via Python AST [REQ-TRC-004]."""

    def test_extract_pytest_markers_and_ignore_helpers(self) -> None:
        py_content = """
import pytest

@pytest.mark.req("REQ-ARCH-006")
def test_zero_copy_numpy_buffer():
    pass

def _helper_not_a_test():
    \"\"\"Helper mentioning REQ-ARCH-001 should not be recorded as a test.\"\"\"
    pass

class TestDiracSuite:
    \"\"\"Suite verifying [REQ-ARCH-001].\"\"\"

    def __init__(self):
        pass

    def _init_fixture(self):
        pass

    def test_hermiticity_check(self):
        pass

    # [REQ-ARCH-003]
    def test_spectral_symmetry(self):
        pass
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            py_file = tmp_root / "test_physics.py"
            py_file.write_text(py_content, encoding="utf-8")

            mappings = vt.extract_python_tests(py_file, tmp_root)
            symbols = {m.full_symbol for m in mappings}

            self.assertIn("test_physics.test_zero_copy_numpy_buffer", symbols)
            self.assertIn("test_physics.TestDiracSuite.test_hermiticity_check", symbols)
            self.assertIn("test_physics.TestDiracSuite.test_spectral_symmetry", symbols)
            # Non-test helpers must be excluded
            self.assertNotIn("test_physics._helper_not_a_test", symbols)
            self.assertNotIn("test_physics.TestDiracSuite.__init__", symbols)
            self.assertNotIn("test_physics.TestDiracSuite._init_fixture", symbols)


class TestJUnitXMLAndPRSummary(unittest.TestCase):
    """Test JUnit XML parsing and PR test-diff detection [REQ-TRC-005]."""

    def test_parse_junit_states(self) -> None:
        xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<testsuites>
  <testsuite name="TestSuite" tests="4">
    <testcase classname="DiracTests" name="PassTest" time="0.001" />
    <testcase classname="DiracTests" name="FailTest" time="0.002">
      <failure message="Expected equality">Mismatch at (0,0)</failure>
    </testcase>
    <testcase classname="DiracTests" name="ErrorTest" time="0.003">
      <error message="std::logic_error">Index out of bounds</error>
    </testcase>
    <testcase classname="DiracTests" name="SkipTest" time="0.000">
      <skipped message="Disabled on macOS" />
    </testcase>
  </testsuite>
</testsuites>
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            xml_file = tmp_root / "junit.xml"
            xml_file.write_text(xml_content, encoding="utf-8")

            results = vt.parse_junit_xml(xml_file)
            self.assertEqual(results["DiracTests.PassTest"].status, "PASSED")
            self.assertEqual(results["DiracTests.FailTest"].status, "FAILED")
            self.assertEqual(results["DiracTests.ErrorTest"].status, "ERROR")
            self.assertEqual(results["DiracTests.SkipTest"].status, "SKIPPED")

    def test_is_test_file_path_detection(self) -> None:
        self.assertTrue(vt.is_test_file_path("src/core/tests/tDiracOperator.cpp"))
        self.assertTrue(vt.is_test_file_path("src/python_bindings/tests/unit/test_dirac.py"))
        self.assertTrue(vt.is_test_file_path("scripts/tests/test_verify_traceability.py"))
        self.assertFalse(vt.is_test_file_path("src/core/DiracOperator.cpp"))
        self.assertFalse(vt.is_test_file_path("docs/eps/ep-5.md"))


class TestShiftLeftGitHooks(unittest.TestCase):
    """Test shift-left local git hooks presence and executable configuration [REQ-TRC-006]."""

    def test_git_hooks_exist_and_are_executable(self) -> None:
        pre_commit = REPO_ROOT / ".githooks" / "pre-commit"
        pre_push = REPO_ROOT / ".githooks" / "pre-push"

        self.assertTrue(pre_commit.exists(), "Missing .githooks/pre-commit")
        self.assertTrue(pre_push.exists(), "Missing .githooks/pre-push")
        self.assertTrue(os.access(pre_commit, os.X_OK), ".githooks/pre-commit must be executable")
        self.assertTrue(os.access(pre_push, os.X_OK), ".githooks/pre-push must be executable")

        pre_commit_text = pre_commit.read_text(encoding="utf-8")
        self.assertIn("verify_traceability.py --audit", pre_commit_text)


class TestCIPlannerAndDependencyGraph(unittest.TestCase):
    """Test C++ #include graph, pybind11 AST bridge, and Two-Hop requirement impact planner [REQ-TRC-005]."""

    def test_pybind11_symbol_map_extraction(self) -> None:
        import ci_plan

        symbol_map = ci_plan.parse_pybind11_symbol_map(REPO_ROOT)
        self.assertEqual(symbol_map.get("DiracOperator"), "src/core/DiracOperator.hpp")
        self.assertEqual(symbol_map.get("Action"), "src/core/BarrettGlaser/Action.hpp")
        self.assertEqual(symbol_map.get("Metropolis"), "src/core/BarrettGlaser/Metropolis.hpp")
        self.assertEqual(symbol_map.get("StdRng"), "src/core/StdRng.hpp")
        self.assertEqual(symbol_map.get("GslRng"), "src/core/StdRng.hpp")
        self.assertEqual(symbol_map.get("set_max_clifford_mode"), "src/core/Clifford.hpp")

    def test_cpp_component_isolation_hamiltonian(self) -> None:
        import ci_plan

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
        import ci_plan

        # Touching Clifford.cpp (which contains zero REQ- tags) must infer REQ-ARCH-001
        # and REQ-ARCH-003 via its dependent C++ unit tests (tClifford.cpp, tDiracOperator.cpp).
        plan = ci_plan.compute_plan_from_files(
            changed_files=["src/core/Clifford.cpp"],
            base_ref="origin/main",
            repo_root=REPO_ROOT,
        )
        self.assertIn("CliffordTests", plan.affected_cpp_suites)
        self.assertIn("DiracOperatorTests", plan.affected_cpp_suites)
        self.assertIn("REQ-ARCH-001", plan.impacted_requirements)
        self.assertIn("REQ-ARCH-003", plan.impacted_requirements)
        self.assertIn("component-test", plan.requirement_triggers["REQ-ARCH-001"])
        # Python tests using rfl.DiracOperator / rfl.set_max_clifford_mode are also affected via pybind11 bridge
        self.assertIn(
            "src/python_bindings/tests/integration/test_boundary_safety.py",
            plan.affected_python_test_files,
        )

    def test_ep_spec_diff_triggers_mapped_tests_without_code_change(self) -> None:
        import ci_plan

        plan = ci_plan.compute_plan_from_files(
            changed_files=["docs/eps/ep-1-core-architecture-redesign.md"],
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
        import ci_plan

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
        import ci_plan

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

    def test_write_github_output(self) -> None:
        import os
        import tempfile
        import ci_plan

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


if __name__ == "__main__":
    unittest.main()
