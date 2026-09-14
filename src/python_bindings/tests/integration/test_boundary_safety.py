"""Integration tests for boundary safety and exception propagation."""
import numpy as np
import pytest
import rfl


def test_exception_bubbling():
    """Verify that invalid parameters raise C++ exceptions in Python."""
    # 1. Verify that a negative matrix dimension raises an exception.
    with pytest.raises(Exception, match="Matrix dimension"):
        rfl.DiracOperator(1, 3, -10)

    # 2. Verify that an oversized Clifford module signature raises an error.
    with pytest.raises(Exception, match="Clifford algebra mode"):
        rfl.DiracOperator(10, 10, 10)

    # 3. Verify that the user can configure the maximum Clifford mode.
    default_max = rfl.get_max_clifford_mode()
    assert default_max == 16

    # Lower the limit to 4 and verify that p+q=5 raises an exception.
    rfl.set_max_clifford_mode(4)
    with pytest.raises(Exception, match="Clifford algebra mode"):
        rfl.DiracOperator(2, 3, 10)

    # Verify that valid parameters succeed under the new limit.
    valid_dirac = rfl.DiracOperator(1, 3, 10)
    assert valid_dirac.get_type() == (1, 3)

    # Restore the default limit.
    rfl.set_max_clifford_mode(default_max)


def test_numpy_integration():
    """Verify that NumPy arrays converted from C++ maintain memory safety."""
    dirac = rfl.DiracOperator(1, 3, 10)
    eigenvals = dirac.get_eigenvalues()

    # Execute standard NumPy operations.
    mean_val = np.mean(eigenvals)
    std_val = np.std(eigenvals)

    # Verify array validity and numeric operations.
    assert isinstance(mean_val, (float, np.floating))
    assert isinstance(std_val, (float, np.floating))
    assert eigenvals.flags.owndata or eigenvals.base is not None


def test_lifecycle_and_gc():
    """Verify object destruction during repeated allocation cycles."""
    for _ in range(1000):
        temp_dirac = rfl.DiracOperator(1, 1, 5)
        temp_eigen = temp_dirac.get_eigenvalues()
        del temp_eigen
        del temp_dirac
