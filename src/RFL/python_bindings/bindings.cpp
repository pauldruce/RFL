#include <armadillo>
#include <carma>
#include <pybind11/pybind11.h>

#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "Clifford.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"

namespace py = pybind11;

PYBIND11_MODULE(rfl, m) {
  m.doc() = "Python bindings for the Random Fuzzy Library (RFL)";

  m.def("set_max_clifford_mode", &Clifford::setMaxMode, py::arg("max_mode"), "Set the maximum allowed Clifford algebra mode (p+q) to prevent excessive memory allocation.");
  m.def("get_max_clifford_mode", &Clifford::getMaxMode, "Get the maximum allowed Clifford algebra mode (p+q)");

  py::class_<IDiracOperator>(m, "IDiracOperator");

  py::class_<DiracOperator, IDiracOperator>(m, "DiracOperator")
      .def(py::init<int, int, int>(), py::arg("p"), py::arg("q"), py::arg("dim"))
      .def("get_type", &DiracOperator::getType)
      .def("get_matrix_dimension", &DiracOperator::getMatrixDimension)
      .def("get_eigenvalues", &DiracOperator::getEigenvalues,
           "Returns the eigenvalues of the Dirac Operator");

  py::class_<Action>(m, "Action")
      .def(py::init<double, double>(), py::arg("g_2"), py::arg("g_4"))
      .def("get_g2", &Action::getG2)
      .def("get_g4", &Action::getG4)
      .def("set_params", &Action::setParams)
      .def("calculate_s", &Action::calculateS);

  py::class_<GslRng>(m, "GslRng")
      .def(py::init<unsigned long>(), py::arg("seed"));

  py::class_<Metropolis>(m, "Metropolis")
      .def(py::init([](double g_2, double g_4, double scale, int num_steps, unsigned long seed) {
             auto action = std::make_unique<Action>(g_2, g_4);
             auto rng = std::make_unique<GslRng>(seed);
             return std::make_unique<Metropolis>(std::move(action), scale, num_steps, std::move(rng));
           }),
           py::arg("g_2"), py::arg("g_4"), py::arg("scale"), py::arg("num_steps"), py::arg("seed"))
      .def("update_dirac", &Metropolis::updateDirac);
}
