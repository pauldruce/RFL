#!/usr/bin/env python3
"""scripts/tests/test_verify_traceability.py

Comprehensive unit test suite for the RFL Requirements Traceability and Governance Engine.
Verifies requirement ID validation, EP table parsing, C++ (GTest/Catch2) extraction,
Python AST extraction, JUnit XML evaluation, manifest hashing, and zero-copy invariants.

Conforms to EP-5 (REQ-TRC-003).
Uses standard Python library only (runnable via unittest or pytest).
"""

from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

# Import the module under test
import sys
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import verify_traceability as vt


class TestRequirementIDValidation(unittest.TestCase):
    """Test requirement ID regex and namespacing rules (REQ-TRC-001)."""

    def test_valid_canonical_ids(self) -> None:
        valid_ids = [
            "REQ-ARCH-001",
            "REQ-PKG-009",
            "REQ-CI-010",
            "REQ-PHY-001",
            "REQ-TRC-008",
            "REQ-AI-042",
        ]
        for req_id in valid_ids:
            self.assertIsNotNone(
                vt.REQ_ID_PATTERN.match(req_id),
                f"Expected valid ID '{req_id}' to match pattern",
            )

    def test_invalid_canonical_ids(self) -> None:
        invalid_ids = [
            "REQ-001",         # Un-namespaced legacy format
            "REQ-A-001",       # Subsystem code too short (< 2 letters)
            "REQ-INVALID-001", # Subsystem code too long (> 4 letters)
            "REQ-arch-001",    # Lowercase subsystem code
            "REQ-ARCH-1",      # Non-zero-padded number
            "REQ-ARCH-01",     # 2-digit number
            "ARCH-001",        # Missing REQ prefix
            "REQ-008b",        # Alphanumeric suffix
        ]
        for req_id in invalid_ids:
            self.assertIsNone(
                vt.REQ_ID_PATTERN.match(req_id),
                f"Expected invalid ID '{req_id}' to fail pattern match",
            )


class TestLegacyRequirementAliasing(unittest.TestCase):
    """Test legacy alias resolution and collision disambiguation (Section 7 Phase 1)."""

    def test_path_prefixed_disambiguation(self) -> None:
        # Disambiguate colliding REQ-001 by path prefix
        arch_id = vt.canonicalize_req("REQ-001", "src/core/tests/tDiracOperator.cpp")
        self.assertEqual(arch_id, "REQ-ARCH-001")

        phy_id = vt.canonicalize_req("REQ-001", "tests/physics/test_invariants.py")
        self.assertEqual(phy_id, "REQ-PHY-001")

    def test_unambiguous_legacy_aliases(self) -> None:
        pkg_id = vt.canonicalize_req("REQ-PKG-01", "any/arbitrary/path.py")
        self.assertEqual(pkg_id, "REQ-PKG-001")

        ci_id = vt.canonicalize_req("REQ-ARC-001", "")
        self.assertEqual(ci_id, "REQ-CI-001")

    def test_already_canonical_ids_unchanged(self) -> None:
        modern_id = vt.canonicalize_req("REQ-ARCH-001", "src/core/tests/tDiracOperator.cpp")
        self.assertEqual(modern_id, "REQ-ARCH-001")


class TestMathPipeShielding(unittest.TestCase):
    """Test LaTeX math pipe character shielding to preserve Markdown table integrity."""

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


class TestEPTableParser(unittest.TestCase):
    """Test Markdown requirement table extraction from Enhancement Proposals."""

    def test_parse_four_column_ep_table(self) -> None:
        ep_content = """# EP-99: Test Proposal

## 3. Research Workflows & Scientific Requirements

### 3.2 Functional Requirements & Invariants

| Requirement ID | Requirement Summary | Lifecycle Status | Operational & Physical Invariant |
| :--- | :--- | :--- | :--- |
| **REQ-ARCH-001** | **Hermiticity Preservation** | Phase 1 (Implemented) | Matrix variations $\\delta M_k$ must preserve $H_i^\\dagger = H_i$. |
| **REQ-ARCH-002** | **Analytic Variation** | Phase 2 (Planned) | Fast $\\mathcal{O}(N^3)$ local variation equals global recomputation. |
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            ep_file = tmp_root / "ep-99-test.md"
            ep_file.write_text(ep_content, encoding="utf-8")

            parsed_requirements, errors = vt.parse_ep_file(ep_file, tmp_root)
            self.assertEqual(len(errors), 0, f"Errors: {errors}")
            self.assertEqual(len(parsed_requirements), 2)

            self.assertEqual(parsed_requirements[0].req_id, "REQ-ARCH-001")
            self.assertEqual(parsed_requirements[0].summary, "Hermiticity Preservation")
            self.assertEqual(parsed_requirements[0].lifecycle_status, "Phase 1 (Implemented)")
            self.assertTrue(parsed_requirements[0].is_implemented)

            self.assertEqual(parsed_requirements[1].req_id, "REQ-ARCH-002")
            self.assertFalse(parsed_requirements[1].is_implemented)


class TestCppTestExtraction(unittest.TestCase):
    """Test requirement tag extraction from GoogleTest and Catch2 source files."""

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
    """Test requirement tag extraction from Pytest files via Python AST."""

    def test_extract_pytest_markers_and_docstrings(self) -> None:
        py_content = """
import pytest

@pytest.mark.req("REQ-ARCH-006")
def test_zero_copy_numpy_buffer():
    pass

def test_hermiticity_check():
    \"\"\"Verify Dirac operator hermiticity [REQ-ARCH-001].\"\"\"
    pass

# [REQ-ARCH-003]
def test_spectral_symmetry():
    pass
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            py_file = tmp_root / "test_physics.py"
            py_file.write_text(py_content, encoding="utf-8")

            mappings = vt.extract_python_tests(py_file, tmp_root)
            self.assertEqual(len(mappings), 3)

            mapped_req_ids = {m.req_id for m in mappings}
            self.assertIn("REQ-ARCH-006", mapped_req_ids)
            self.assertIn("REQ-ARCH-001", mapped_req_ids)
            self.assertIn("REQ-ARCH-003", mapped_req_ids)


class TestManifestVerification(unittest.TestCase):
    """Test cryptographic SHA-256 integrity verification of sealed qualification tests."""

    def test_manifest_verification_normalises_line_endings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            sealed_file = tmp_root / "test_sealed.py"
            content = "def test_sealed_physics():\n    assert True\n"
            sealed_file.write_bytes(content.encode("utf-8"))

            expected_hash = hashlib.sha256(content.encode("utf-8")).hexdigest()
            manifest_file = tmp_root / ".manifest"
            manifest_file.write_text(f"{expected_hash}  test_sealed.py\n", encoding="utf-8")

            success, errors = vt.verify_manifest(manifest_file, tmp_root)
            self.assertTrue(success, f"Manifest verification failed: {errors}")

    def test_manifest_verification_detects_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            sealed_file = tmp_root / "test_sealed.py"
            sealed_file.write_bytes(b"original content\n")

            manifest_file = tmp_root / ".manifest"
            manifest_file.write_text(f"{'0' * 64}  test_sealed.py\n", encoding="utf-8")

            success, errors = vt.verify_manifest(manifest_file, tmp_root)
            self.assertFalse(success)
            self.assertIn("Integrity violation", errors[0])


class TestJUnitXMLParser(unittest.TestCase):
    """Test JUnit XML parsing across PASSED, FAILED, ERROR, and SKIPPED states."""

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


class TestZeroCopyFortranInvariant(unittest.TestCase):
    """Test REQ-TRC-008 Fortran contiguous buffer verification helper."""

    def test_fortran_contiguous_verification(self) -> None:
        class MockArray:
            def __init__(self, f_contiguous: bool, ptr: int = 0x1000):
                self.flags = {"F_CONTIGUOUS": f_contiguous}
                class CTypes:
                    pass
                self.ctypes = CTypes()
                self.ctypes.data = ptr

        valid_arr = MockArray(f_contiguous=True, ptr=0x2000)
        self.assertTrue(vt.assert_zero_copy_fortran_aligned(valid_arr, 0x2000))

        invalid_arr = MockArray(f_contiguous=False)
        with self.assertRaises(AssertionError):
            vt.assert_zero_copy_fortran_aligned(invalid_arr)

        with self.assertRaises(AssertionError):
            vt.assert_zero_copy_fortran_aligned(valid_arr, 0x9999)


if __name__ == "__main__":
    unittest.main()
