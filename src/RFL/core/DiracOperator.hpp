//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_DIRACOPERATOR_HPP
#define RFL_DIRACOPERATOR_HPP

#include "Clifford.hpp"
#include "IDiracOperator.hpp"
#include "IRng.hpp"
#include <armadillo>
#include <memory>

/**
 * @class DiracOperator
 *
 * @brief Represents a Dirac operator with specified signature and matrix dimension.
 *
 * Provides methods to access properties and calculate quantities related to the Dirac operator.
 */
class DiracOperator final : public IDiracOperator {
public:
  /**
   * Constructs a Dirac operator with signature (p,q) and matrix dimension dim.
   *
   * @param p Number of Hermitian gamma matrices.
   * @param q Number of anti-Hermitian gamma matrices.
   * @param dim Matrix dimension of the H and L matrices.
   */
  DiracOperator(int p, int q, int dim);

  /**
   * Copy constructor for DiracOperator.
   *
   * @param original DiracOperator object to copy.
   */
  DiracOperator(const DiracOperator& original);

  /**
   * Returns the signature of the Dirac operator as a pair (p,q).
   *
   * @return std::pair<int,int> containing (p,q).
   */
  std::pair<int, int> getType() const override { return std::pair{m_clifford.getP(), m_clifford.getQ()}; }

  /**
   * Returns the matrix dimension of the H and L matrices.
   */
  int getMatrixDimension() const override { return m_dim; };

  /**
   * Returns the dimension of the gamma matrices.
   */
  int getGammaDimension() const override { return m_gamma_dim; };

  /**
   * Returns the total number of H and L matrices.
   */
  int getNumMatrices() const override { return m_num_matrices; };

  /**
   * Returns the number of Hermitian matrices (H).
   */
  int getNumHermitianMatrices() const override { return m_num_herm; };

  /**
   * Returns the number of anti-Hermitian matrices (L).
   */
  int getNumAntiHermitianMatrices() const override { return m_num_antiherm; };

  /**
   * Returns a reference to the vector of H and L matrices.
   */
  std::vector<arma::cx_mat>& getMatrices() const override { return *m_matrices; }

  /**
   * Returns a reference to the vector of signs (+1 or -1).
   *
   * The sign matches the matrix at the same index in getMatrices().
   * A value of +1 indicates a Hermitian matrix (H).
   * A value of -1 indicates an anti-Hermitian matrix (L).
   */
  std::vector<int>& getEpsilons() const override { return *m_epsilons; }

  /**
   * Returns a reference to the vector of conjugate momenta matrices.
   */
  // TODO: Move to Hamiltonian class or document usage.
  std::vector<arma::cx_mat>& getMomenta() const override { return *m_momenta; }

  /**
   * Calculates the trace of the Dirac operator squared, \f$\text{Tr}(D^2)\f$.
   *
   * @return Trace value as a double.
   */
  double traceOfDiracSquared() const override;

  /**
   * Calculates the trace of the fourth power of the Dirac operator, \f$\text{Tr}(D^4)\f$.
   *
   * @return Trace value as a double.
   */
  double traceOfDirac4() const override;

  /**
   * Constructs and returns the assembled matrix representation of the Dirac operator.
   *
   * @return Assembled Dirac operator matrix.
   */
  arma::cx_mat getDiracMatrix() const override;

  /**
   * Calculates and returns the eigenvalues of the assembled Dirac operator matrix.
   *
   * @return Armadillo vector of real eigenvalues.
   */
  arma::vec getEigenvalues() const override;

  /**
   * Returns a vector of the Hermitian matrices (H_i) in the Dirac operator.
   */
  std::vector<arma::cx_mat> getHermitianMatrices() const override;

  /**
   * Returns a vector of the anti-Hermitian matrices (L_j) in the Dirac operator.
   */
  std::vector<arma::cx_mat> getAntiHermitianMatrices() const override;

  /**
   * Calculates the derivative of the Barrett-Glaser action with respect to matrix element k.
   *
   * @param k Index of the matrix.
   * @param herm True if the matrix is Hermitian, false otherwise.
   * @param g_2 Coupling constant for Tr(D^2).
   * @return Derivative matrix.
   */
  arma::cx_mat derDirac24(const int& k, const bool& herm, double g_2) const;

  /**
   * Calculates the derivative of Tr(D^2) with respect to matrix element k.
   *
   * @param k Index of the matrix.
   * @return Derivative matrix.
   */
  arma::cx_mat derDirac2(const int& k) const;

  /**
   * Calculates the derivative of Tr(D^4) with respect to matrix element k.
   *
   * @param k Index of the matrix.
   * @param herm True if the matrix is Hermitian, false otherwise.
   * @return Derivative matrix.
   */
  arma::cx_mat derDirac4(const int& k, const bool& herm) const;

  /**
   * Returns a reference to the four-product table of omega matrices.
   */
  std::vector<arma::cx_double>& getOmegaTable4() const override { return *m_omega_table_4; }

  /**
   * Prints the four-product table of omega matrices to standard output.
   */
  void printOmegaTable4() const;

  /**
   * Randomises all internal matrix elements using the provided random number generator.
   *
   * @param rng Random number generator engine.
   */
  void randomiseMatrices(const IRng& rng) const override;

private:
  // Clifford module for the Dirac operator.
  Clifford m_clifford;
  // Matrix dimension of the H and L matrices.
  int m_dim;
  // Total number of matrices, number of Hermitian matrices, and number of anti-Hermitian matrices.
  int m_num_matrices, m_num_herm, m_num_antiherm;
  // Dimension of gamma matrices.
  int m_gamma_dim;
  // Internal H and L matrices (all Hermitian representations).
  std::unique_ptr<std::vector<arma::cx_mat>> m_matrices;
  // Conjugate momenta matrices.
  std::unique_ptr<std::vector<arma::cx_mat>> m_momenta;
  // Omega matrices.
  std::unique_ptr<std::vector<arma::cx_mat>> m_omegas;
  // Sign vector: +1 for Hermitian (H), -1 for anti-Hermitian (L).
  std::unique_ptr<std::vector<int>> m_epsilons;
  // Four-product table for omega matrices.
  std::unique_ptr<std::vector<arma::cx_double>> m_omega_table_4{};

  void initOmegaTable4();
  double computeA4(const int& i_1, const int& i_2, const int& i_3, const int& i_4) const;
  double computeA2(const int& i_1, const int& i_2) const;
  double computeA(const int& i) const;
  arma::cx_mat computeB4(const int& k,
                         const int& i_2,
                         const int& i_3,
                         const int& i_4,
                         const double& cliff,
                         const bool& neg) const;
  arma::cx_mat computeB2(const int& k, const int& i) const;
  arma::cx_mat computeB(const int& k) const;
};
#endif// RFL_DIRACOPERATOR_HPP
