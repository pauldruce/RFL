//
// Created by Paul Druce on 30/12/2023.
//

#ifndef EXAMPLES_EIGENVALUERECORDER_HPP
#define EXAMPLES_EIGENVALUERECORDER_HPP

#define ARMA_USE_HDF5
#include "DiracOperator.hpp"
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <utility>

class EigenvalueRecorder {
public:
  EigenvalueRecorder(const DiracOperator& dirac, double g2, std::string simulationId)
      : m_dirac(dirac), m_g2(g2), m_simulationId(std::move(simulationId)) {
  }

  void recordEigenvalues(int diracId);

private:
  const DiracOperator& m_dirac;
  const double m_g2;
  std::string m_simulationId;
};

#endif//EXAMPLES_EIGENVALUERECORDER_HPP