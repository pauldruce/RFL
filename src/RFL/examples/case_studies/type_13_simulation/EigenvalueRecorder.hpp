//
// Created by Paul Druce on 30/12/2023.
//

#ifndef EXAMPLES_EIGENVALUERECORDER_HPP
#define EXAMPLES_EIGENVALUERECORDER_HPP

#define ARMA_USE_HDF5
#include "DiracOperator.hpp"

/**
 * @class EigenvalueRecorder
 * @brief Records the eigenvalue spectrum of a Dirac operator to an HDF5 file.
 *
 * This recorder uses an HDF5 file to store eigenvalues.
 *
 * The environment variable `RFL_OUTPUT_DIR` controls the output directory.
 * If this variable is not set, the recorder writes output to `/tmp/RFL`.
 *
 * Each recorder writes to an HDF5 dataset group named by `simulationId`.
 * The recorder appends eigenvalue datasets to this group.
 */
class EigenvalueRecorder {
public:
  EigenvalueRecorder(const DiracOperator& dirac, const double g2, const std::string& simulationId)
      : m_dirac(dirac), m_g2(g2), m_simulationId(simulationId) {
  }

  /**
   * @brief Records the eigenvalue spectrum of the Dirac operator.
   *
   * @param diracId Identifier for the Dirac operator configuration.
   */
  void recordEigenvalues(int diracId) const;

private:
  const DiracOperator& m_dirac;
  const double m_g2;
  std::string m_simulationId;
};

#endif//EXAMPLES_EIGENVALUERECORDER_HPP