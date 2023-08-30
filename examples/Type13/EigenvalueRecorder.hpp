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

/**
 * @class EigenvalueRecorder
 * @brief A simple class to record the eigenvalues of a DiracOperator to a file.
 *
 * This recorder uses a HDF5 file to store the eigenvalues. The file name/location is
 * set by the parameters of the simulation - specifically the type of the DiracOperator and
 * the g2 value of action.
 * Each recorder writes to a dataset group in the HDF5 file named by the constructor
 * parameter "simulationId". Each eigenvalue data set is appended to this dataset group.
 *
 */
class EigenvalueRecorder {
public:
  EigenvalueRecorder(const DiracOperator& dirac, double g2, std::string simulationId)
      : m_dirac(dirac), m_g2(g2), m_simulationId(std::move(simulationId)) {
  }

  /**
   * @brief Records the eigenvalues of a given Dirac operator.
   *
   * The recordEigenvalues function accepts a diracId as input and records the eigenvalues of
   * the dirac operator to a file.
   */
  void recordEigenvalues(int diracId);

private:
  const DiracOperator& m_dirac;
  const double m_g2;
  std::string m_simulationId;
};

#endif//EXAMPLES_EIGENVALUERECORDER_HPP