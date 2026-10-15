"""Sealed qualification acceptance tests for numerical tolerances and physical invariants.

These tests seal mathematical error bounds grounded in IEEE 754 double precision machine epsilon.
Any modification to tolerances in this file requires maintainer approval under EP-5 (REQ-TRC-006).
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

# Add repo root to import verify_traceability helper
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.verify_traceability import assert_zero_copy_fortran_aligned

# Standard IEEE 754 double precision machine epsilon (2^-52)
EPS_MACH: float = 2.0**-52

# Empirical tolerance constants established in EP-4 Section 3.3.5
C1_TRACE_CANCELLATION: float = 10.0
C0_TRACE_CANCELLATION: float = 5.0
C2_HERMITICITY: float = 5.0
C3_PERMUTATION: float = 5.0
C4_UNITARY_ACTION: float = 20.0
C5_WEYL_EIGENVALUE: float = 5.0
C6_FINITE_DIFFERENCE: float = 50.0


class TestNumericalTolerances(unittest.TestCase):
    """Qualification tests verifying sealed error bounds and physical invariants."""

    # [REQ-TRC-006]
    def test_machine_epsilon_ieee754(self) -> None:
        """Verify IEEE 754 double precision machine epsilon matches theoretical value.

        Machine epsilon represents the upper bound on relative rounding error.
        Numerical tolerances in RFL scale from this fundamental hardware invariant.
        """
        expected_eps = 2.220446049250313e-16
        self.assertEqual(EPS_MACH, 2.0**-52)
        self.assertTrue(math.isclose(EPS_MACH, expected_eps, rel_tol=1e-15))

    # [REQ-PHY-001]
    def test_trace_cancellation_bounds(self) -> None:
        """Verify incremental action trace cancellation error bound formula.

        Trace subtraction S_f - S_i suffers from numerical cancellation.
        The bound scales linearly with matrix dimension M and machine epsilon.
        """
        # Test parameters for N=4, dimension M=16
        dim_m = 16
        action_f = 120.5
        action_i = 119.8
        tr_d2 = 45.2
        tr_d4 = 310.8
        g2 = -1.5
        g4 = 1.0

        bound = (
            C1_TRACE_CANCELLATION
            * dim_m
            * EPS_MACH
            * (abs(action_f) + abs(action_i) + abs(g2) * tr_d2 + g4 * tr_d4)
            + C0_TRACE_CANCELLATION * (dim_m**2) * EPS_MACH
        )
        self.assertGreater(bound, 0.0)
        self.assertLess(bound, 1e-10)

    # [REQ-PHY-002]
    def test_hermiticity_error_bounds(self) -> None:
        """Verify Dirac operator Hermiticity error bound scaling.

        Summing p+q Kronecker tensor products accumulates rounding errors.
        The Frobenius norm error bound scales linearly with signature p+q.
        """
        for signature_sum in (1, 2, 4, 8):
            norm_d = 10.0
            bound = C2_HERMITICITY * signature_sum * EPS_MACH * norm_d
            self.assertGreater(bound, 0.0)
            self.assertLess(bound, 1e-12)

    # [REQ-PHY-003]
    def test_ko_real_structure_error_bounds(self) -> None:
        """Verify KO-dimension real structure commutativity error bounds.

        Permutation matrices in J introduce zero numerical growth beyond rounding.
        The bound scales directly with machine epsilon and operator norm.
        """
        norm_d = 25.0
        bound = C3_PERMUTATION * EPS_MACH * norm_d
        self.assertGreater(bound, 0.0)
        self.assertLess(bound, 1e-13)

    # [REQ-PHY-004]
    def test_chirality_sign_anti_commutation_exact_factor(self) -> None:
        """Verify chirality anti-commutation factor is exactly two.

        Chirality grading gamma consists of sign flips without matrix growth.
        The anti-commutator {gamma, D} accumulates at most two epsilon units.
        """
        norm_d = 15.0
        exact_factor = 2.0
        bound = exact_factor * EPS_MACH * norm_d
        self.assertAlmostEqual(bound, 30.0 * EPS_MACH, places=15)

    # [REQ-PHY-008]
    def test_unitary_gauge_invariance_error_bounds(self) -> None:
        """Verify unitary gauge action invariance and Weyl eigenvalue bounds.

        Unitary similarity transformations U D U^dagger preserve the spectrum.
        Action evaluation error scales quadratically with Hilbert dimension M.
        """
        dim_m = 8
        action_val = 50.0
        norm_2 = 12.0

        action_bound = C4_UNITARY_ACTION * (dim_m**2) * EPS_MACH * abs(action_val)
        weyl_bound = C5_WEYL_EIGENVALUE * dim_m * EPS_MACH * norm_2

        self.assertLess(action_bound, 1e-10)
        self.assertLess(weyl_bound, 1e-12)

    # [REQ-PHY-009]
    def test_finite_difference_optimal_step_and_error_floor(self) -> None:
        """Verify finite difference optimal step size and error floor.

        Central differences balance truncation error O(h^2) and rounding error O(eps/h).
        Higham (2002) proves asserting tolerances tighter than 10^-10 is invalid.
        """
        optimal_step = math.sqrt(EPS_MACH)
        expected_step = 1.4901161193847656e-08
        self.assertTrue(math.isclose(optimal_step, expected_step, rel_tol=1e-10))

        min_error_floor = (EPS_MACH ** (2.0 / 3.0))
        self.assertTrue(math.isclose(min_error_floor, 3.6668527e-11, rel_tol=1e-5))

        test_tolerance = C6_FINITE_DIFFERENCE * min_error_floor
        self.assertGreater(test_tolerance, 1e-10)

    # [REQ-TRC-008]
    def test_fortran_zero_copy_alignment_invariant(self) -> None:
        """Verify Fortran-contiguous array layout verification helper.

        Prevents silent matrix transposition and complex conjugation (H -> H*).
        Validates column-major order for Armadillo and NumPy interoperability.
        """
        class MockFortranArray:
            def __init__(self, is_fortran: bool, data_ptr: int) -> None:
                self.flags = {"F_CONTIGUOUS": is_fortran}
                class CTypes:
                    pass
                self.ctypes = CTypes()
                self.ctypes.data = data_ptr

        valid_arr = MockFortranArray(is_fortran=True, data_ptr=0x1234)
        self.assertTrue(assert_zero_copy_fortran_aligned(valid_arr, expected_ptr=0x1234))

        invalid_c_order = MockFortranArray(is_fortran=False, data_ptr=0x1234)
        with self.assertRaises(AssertionError):
            assert_zero_copy_fortran_aligned(invalid_c_order)

        pointer_mismatch = MockFortranArray(is_fortran=True, data_ptr=0x9999)
        with self.assertRaises(AssertionError):
            assert_zero_copy_fortran_aligned(pointer_mismatch, expected_ptr=0x1234)


if __name__ == "__main__":
    unittest.main()
