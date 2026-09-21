"""Random Fuzzy Library (RFL) Python package.

Provides high-performance Markov Chain Monte Carlo simulations of Finite
Noncommutative Geometries. GslRng is deprecated and maintained as an alias
for StdRng.
"""

from ._rfl import *  # noqa: F403
from ._rfl import __version__

__all__ = [
    "Action",
    "DiracOperator",
    "GslRng",
    "IDiracOperator",
    "Metropolis",
    "StdRng",
    "__version__",
    "get_max_clifford_mode",
    "set_max_clifford_mode",
]
