//
// Created by Paul Druce on 30/12/2023.
//

#ifndef EXAMPLES_EIGENVALUERECORDER_HPP
#define EXAMPLES_EIGENVALUERECORDER_HPP

#define ARMA_USE_HDF5
#include "DiracOperator.hpp"

/**
 * @class EigenvalueRecorder
 * @brief A simple class to record the eigenvalues of a DiracOperator to a file.
 *
 * This recorder uses a HDF5 file to store the eigenvalues.
 *
 * The output directory of the files is controlled by the environment variable
 * RFL_OUTPUT_DIR.
 * If this variable is not set, then all output is placed in /tmp/RFL
 *
 * Each recorder writes to a dataset group in the HDF5 file named by the constructor
 * parameter "simulationId". Each eigenvalue data set is appended to this dataset group.
 *
 */
class EigenvalueRecorder {
public:
  EigenvalueRecorder(const DiracOperator& dirac, const double g2, const std::string& simulationId)
      : m_dirac(dirac), m_g2(g2), m_simulationId(simulationId) {
  }

  /**
   * @brief Records the eigenvalues of a given Dirac operator.
   *
   * The recordEigenvalues function accepts a diracId as input and records the eigenvalues of
   * the dirac operator to a file.
   */
  void recordEigenvalues(int diracId) const;

private:
  const DiracOperator& m_dirac;
  const double m_g2;
  std::string m_simulationId;
};

#endif//EXAMPLES_EIGENVALUERECORDER_HPP