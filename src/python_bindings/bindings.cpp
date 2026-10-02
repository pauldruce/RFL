// clang-format off
#ifndef CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
#define CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
#endif
#include <carma>
#include <armadillo>
// clang-format on
#include <pybind11/pybind11.h>

#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "Clifford.hpp"
#include "DiracOperator.hpp"
#include "StdRng.hpp"

#ifndef RFL_VERSION_STRING
#define RFL_VERSION_STRING "0.0.0-devel"
#endif

namespace py = pybind11;

PYBIND11_MODULE(_rfl, m) {
  m.doc() = "Python bindings for the Random Fuzzy Library (RFL).";
  m.attr("__version__") = RFL_VERSION_STRING;

  m.def("set_max_clifford_mode", &Clifford::setMaxMode, py::arg("max_mode"), "Set the maximum allowed Clifford algebra mode (p+q) to prevent excessive memory allocation.");
  m.def("get_max_clifford_mode", &Clifford::getMaxMode, "Get the maximum allowed Clifford algebra mode (p+q).");

  py::class_<IDiracOperator>(m, "IDiracOperator");

  py::class_<DiracOperator, IDiracOperator>(m, "DiracOperator")
      .def(py::init<int, int, int>(), py::arg("p"), py::arg("q"), py::arg("dim"))
      .def("get_type", &DiracOperator::getType)
      .def("get_matrix_dimension", &DiracOperator::getMatrixDimension)
      .def("get_eigenvalues", &DiracOperator::getEigenvalues,
           "Return the eigenvalue spectrum of the Dirac operator.");

  py::class_<Action>(m, "Action",
                     R"doc(
Barrett-Glaser spectral action evaluator.

Evaluates the spectral action functional:
    S(D) = g_2 Tr(D^2) + g_4 Tr(D^4)

Parameters
----------
g_2 : float
    Quadratic coupling constant.
g_4 : float
    Quartic coupling constant.

References
----------
Barrett, J. W. & Glaser, L. (2016). arXiv:1510.01377.
See `docs/theory/barrett_glaser.md` for mathematical derivations.
)doc")
      .def(py::init<double, double>(), py::arg("g_2"), py::arg("g_4"),
           "Construct an Action with quadratic (g_2) and quartic (g_4) couplings.")
      .def("get_g2", &Action::getG2, "Return the quadratic coupling constant g_2.")
      .def("get_g4", &Action::getG4, "Return the quartic coupling constant g_4.")
      .def("set_params", &Action::setParams, py::arg("g_2"), py::arg("g_4"),
           "Set quadratic (g_2) and quartic (g_4) coupling constants.")
      .def("calculate_s", &Action::calculateS, py::arg("dirac"),
           R"doc(
Calculate the Barrett-Glaser action using component matrix traces.

Parameters
----------
dirac : IDiracOperator
    The Dirac operator state.

Returns
-------
float
    The calculated spectral action value S(D).

Notes
-----
Evaluates S(D) = g_2 Tr(D^2) + g_4 Tr(D^4) using Clifford trace
factorisations without assembling the full Dirac operator matrix.
)doc")
      .def("calculate_s_from_dirac", &Action::calculateSFromDirac, py::arg("dirac"),
           R"doc(
Calculate the Barrett-Glaser action from the assembled Dirac operator matrix.

Parameters
----------
dirac : IDiracOperator
    The Dirac operator state.

Returns
-------
float
    The calculated spectral action value S(D).

Notes
-----
Assembles the full Dirac matrix D and evaluates S(D) = g_2 Re(Tr(D^2)) + g_4 Re(Tr(D^4)).
Used primarily as a brute-force baseline to verify calculate_s.
)doc");

  // Canonical class name alias
  m.attr("BarrettGlaserAction") = m.attr("Action");

  py::class_<StdRng>(m, "StdRng")
      .def(py::init<>(), "Construct an unseeded StdRng engine.")
      .def(py::init<uint64_t>(), py::arg("seed"), "Construct a seeded StdRng engine.")
      .def("get_uniform", &StdRng::getUniform, "Draw uniform random number in [0, 1).")
      .def("get_gaussian", &StdRng::getGaussian, py::arg("sigma"), "Draw Gaussian random variable.")
      .def("get_uniform_int", &StdRng::getUniformInt, py::arg("min"), py::arg("max"), "Draw discrete uniform integer in [min, max].");

  // Backwards compatibility alias for deprecated GslRng
  m.attr("GslRng") = m.attr("StdRng");

  py::class_<Metropolis>(m, "Metropolis")
      .def(py::init([](double g_2, double g_4, double scale, int num_steps, uint64_t seed) {
             auto action = std::make_unique<Action>(g_2, g_4);
             auto rng = std::make_unique<StdRng>(seed);
             return std::make_unique<Metropolis>(std::move(action), scale, num_steps, std::move(rng));
           }),
           py::arg("g_2"), py::arg("g_4"), py::arg("scale"), py::arg("num_steps"), py::arg("seed"))
      .def("update_dirac", &Metropolis::updateDirac);
}
