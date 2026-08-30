"""Unit tests for the RFL Python core API."""

import numpy as np
import pytest
import rfl


def test_dirac_operator():
    """Verify Dirac operator initialisation, type, and eigenvalue spectrum."""
    p, q = 1, 3
    dim = 10
    dirac = rfl.DiracOperator(p, q, dim)

    assert dirac.get_type() == (p, q)
    assert dirac.get_matrix_dimension() == dim

    eigenvals = dirac.get_eigenvalues()
    assert isinstance(eigenvals, np.ndarray)
    assert len(eigenvals) == 400


def test_action():
    """Verify Barrett-Glaser action parameters and action calculation."""
    g_2 = -1.0
    g_4 = 1.0
    action = rfl.Action(g_2, g_4)

    assert action.get_g2() == g_2
    assert action.get_g4() == g_4

    dirac = rfl.DiracOperator(1, 3, 10)
    s = action.calculate_s(dirac)
    assert isinstance(s, float)


def test_gsl_rng():
    """Verify GSL random number generator initialisation."""
    rng = rfl.GslRng(42)
    assert rng is not None


def test_metropolis():
    """Verify Metropolis algorithm execution and acceptance rate output."""
    dirac = rfl.DiracOperator(1, 3, 10)

    g_2 = -1.0
    g_4 = 1.0
    scale = 1.0
    num_steps = 10
    seed = 42

    metropolis = rfl.Metropolis(g_2, g_4, scale, num_steps, seed)
    acceptance_rate = metropolis.update_dirac(dirac)

    assert 0.0 <= acceptance_rate <= 1.0


