"""Public API surface and signature assurance tests for rfl."""

import packaging.version
import pytest
import rfl


EXPECTED_PUBLIC_SYMBOLS = {
    "Action",
    "DiracOperator",
    "GslRng",
    "IDiracOperator",
    "Metropolis",
    "StdRng",
    "__version__",
    "get_max_clifford_mode",
    "set_max_clifford_mode",
}


def test_public_all_exports():
    """Verify that __all__ contains the exact intended public symbols."""
    assert set(rfl.__all__) == EXPECTED_PUBLIC_SYMBOLS
    assert len(rfl.__all__) == len(EXPECTED_PUBLIC_SYMBOLS)


def test_public_symbols_exist_and_accessible():
    """Verify that each exported symbol is accessible on the module."""
    for symbol in EXPECTED_PUBLIC_SYMBOLS:
        assert hasattr(rfl, symbol), f"Symbol '{symbol}' missing from rfl module."
        attr = getattr(rfl, symbol)
        assert attr is not None


def test_no_third_party_or_internal_leakage():
    """Verify that internal dependencies do not leak into public symbols."""
    disallowed_substrings = ["carma", "numpy", "pybind11", "scikit"]
    for symbol in rfl.__all__:
        for disallowed in disallowed_substrings:
            assert disallowed not in symbol.lower(), f"Leaked internal symbol: {symbol}"

    module_dir = dir(rfl)
    assert "carma" not in module_dir
    assert "numpy" not in module_dir


def test_clifford_mode_functions_signature():
    """Verify Clifford helper function signatures and keyword arguments."""
    assert callable(rfl.set_max_clifford_mode)
    assert callable(rfl.get_max_clifford_mode)

    doc_set = rfl.set_max_clifford_mode.__doc__ or ""
    assert "set_max_clifford_mode(max_mode: int)" in doc_set

    doc_get = rfl.get_max_clifford_mode.__doc__ or ""
    assert "get_max_clifford_mode()" in doc_get

    rfl.set_max_clifford_mode(max_mode=9)
    assert rfl.get_max_clifford_mode() == 9


def test_dirac_operator_surface():
    """Verify DiracOperator constructor signature, methods, and hierarchy."""
    assert issubclass(rfl.DiracOperator, rfl.IDiracOperator)

    doc = rfl.DiracOperator.__init__.__doc__ or ""
    assert "p: int" in doc
    assert "q: int" in doc
    assert "dim: int" in doc

    dirac = rfl.DiracOperator(p=1, q=3, dim=6)
    assert dirac.get_type() == (1, 3)
    assert dirac.get_matrix_dimension() == 6

    eigenvalues = dirac.get_eigenvalues()
    assert hasattr(eigenvalues, "shape")


def test_action_surface():
    """Verify Action constructor signature, methods, and parameters."""
    doc = rfl.Action.__init__.__doc__ or ""
    assert "g_2: float" in doc
    assert "g_4: float" in doc

    action = rfl.Action(g_2=-1.5, g_4=1.0)
    assert action.get_g2() == -1.5
    assert action.get_g4() == 1.0

    dirac = rfl.DiracOperator(p=1, q=3, dim=6)
    value = action.calculate_s(dirac)
    assert isinstance(value, float)


def test_std_rng_surface():
    """Verify StdRng and alias GslRng constructor signatures and sampling methods."""
    doc = rfl.StdRng.__init__.__doc__ or ""
    assert "seed: int" in doc

    # Seeded instantiation
    rng = rfl.StdRng(seed=12345)
    assert rng is not None

    # Unseeded instantiation
    unseeded = rfl.StdRng()
    assert unseeded is not None

    # Sampling method checks
    u = rng.get_uniform()
    assert isinstance(u, float)
    assert 0.0 <= u < 1.0

    g = rng.get_gaussian(sigma=1.0)
    assert isinstance(g, float)

    k = rng.get_uniform_int(min=1, max=10)
    assert isinstance(k, int)
    assert 1 <= k <= 10

    # Backwards compatibility alias check
    assert rfl.GslRng is rfl.StdRng
    compat_rng = rfl.GslRng(seed=42)
    assert compat_rng is not None


def test_metropolis_surface():
    """Verify Metropolis constructor signature and update method."""
    doc = rfl.Metropolis.__init__.__doc__ or ""
    for param in ["g_2: float", "g_4: float", "scale: float", "num_steps: int", "seed: int"]:
        assert param in doc, f"Parameter '{param}' missing from Metropolis.__init__ docstring."

    sampler = rfl.Metropolis(g_2=-1.0, g_4=1.0, scale=0.5, num_steps=2, seed=42)
    dirac = rfl.DiracOperator(p=1, q=3, dim=6)
    rate = sampler.update_dirac(dirac)
    assert 0.0 <= rate <= 1.0


def test_version_string_pep440():
    """Verify package version format conforms to PEP 440."""
    version_str = rfl.__version__
    assert isinstance(version_str, str)
    parsed = packaging.version.parse(version_str)
    assert parsed is not None
    assert parsed.base_version != ""
