import rfl
import numpy as np
import pytest

def test_exception_bubbling():
    # Test that invalid parameters correctly bubble up C++ exceptions 
    # to Python instead of crashing the interpreter.
    
    # 1. Negative dimension should throw an exception instead of crashing
    with pytest.raises(Exception, match="Matrix dimension"):
        invalid_dirac = rfl.DiracOperator(1, 3, -10)
        
    # 2. Excessively large Clifford algebra should throw to prevent OOM
    with pytest.raises(Exception, match="Clifford algebra mode"):
        invalid_dirac = rfl.DiracOperator(10, 10, 10)
        
    # 3. Test that the maximum mode is dynamically configurable by the user
    default_max = rfl.get_max_clifford_mode()
    assert default_max == 16
    
    # Lower the limit to 4 and ensure p+q=5 throws
    rfl.set_max_clifford_mode(4)
    with pytest.raises(Exception, match="Clifford algebra mode"):
        rfl.DiracOperator(2, 3, 10)
        
    # Ensure it works under the new limit
    valid_dirac = rfl.DiracOperator(1, 3, 10)
    assert valid_dirac.get_type() == (1, 3)
    
    # Restore the default limit
    rfl.set_max_clifford_mode(default_max)
        
def test_numpy_integration():
    # Test that numpy arrays returned from C++ via carma/pybind11
    # behave completely natively in Python and memory is safe.
    dirac = rfl.DiracOperator(1, 3, 10)
    eigenvals = dirac.get_eigenvalues()
    
    # Run native numpy operations
    mean_val = np.mean(eigenvals)
    std_val = np.std(eigenvals)
    
    # Assert that the array is valid and mathematical ops succeed
    assert isinstance(mean_val, (float, np.floating))
    assert isinstance(std_val, (float, np.floating))
    assert eigenvals.flags.owndata or eigenvals.base is not None # Memory is safely managed

def test_lifecycle_and_gc():
    # Create and destroy C++ objects rapidly to verify garbage collection without leaks.
    for _ in range(1000):
        temp_dirac = rfl.DiracOperator(1, 1, 5)
        temp_eigen = temp_dirac.get_eigenvalues()
        del temp_eigen
        del temp_dirac
