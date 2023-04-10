//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_DIRACOPERATOR_HPP
#define RFL_DIRACOPERATOR_HPP

#include "Clifford.hpp"
#include "IRng.hpp"
#include <armadillo>
#include <memory>

/**
 * A class that represents the Dirac operator as described in the paper:

 * John W. Barrett, "Matrix geometries and fuzzy spaces as finite spectral triples",
 * J. Math. Phys. 56, 082301 (2015) https://doi.org/10.1063/1.4927224
 *
 */
class DiracOperator {
public:
  /**
   * The default constructor for this class has been deleted, see the other options
   * that are available
   */
  DiracOperator() = delete;
  ~DiracOperator();

  /**
   * The constructor for this class. The parameters p and q represent the
   * Clifford type for the Dirac operator \f$(p,q)\f$. The parameter dim is the
   * dimension of the H and L matrices for the Dirac operator.
   */
  DiracOperator(int p, int q, int dim);

  /**
   * getMatrixDimension returns the dimension of the H and L matrices of the
   * Dirac operator
   */
  int getMatrixDimension() const { return m_dim; };

  /**
   * getGammaDimension returns the dimension of the gamma matrices of the Dirac
   * operator
   */
  int getGammaDimension() const { return m_gamma_dim; };

  /**
   * getNumMatrices returns the total number of H and L matrices
   */
  int getNumMatrices() const { return m_num_matrices; };

  /**
   * getNumHermitianMatrices returns the number of H matrices
   */
  int getNumHermitianMatrices() const { return m_num_herm; };

  /**
   * getNumHermitianMatrices returns the number of L matrices
   */
  int getNumAntiHermitianMatrices() const { return m_num_antiherm; };

  /**
   * This method returns a reference to the vector of H and L matrices of the Dirac operator
   */
  std::vector<arma::cx_mat>& getMatrices() const { return *m_matrices; }

  /**
   * getEpsilons returns a reference to a vector of +/-1. The sign of the entry relates
   * to the matrix with the same index returned from getMatrices. If the value is +1 then
   * the associated matrix is Hermitian. If the value is -1, then the associated matrix
   * is anti-Hermitian.
   */
  std::vector<int>& getEpsilons() const { return *m_epsilons; }

  /**
   * getMomenta returns a reference to a vector of matrices which correspond to..
   */
  // TODO: This should probably be moved to the Hamiltonian class, or needs good documentation about it's use
  std::vector<arma::cx_mat>& getMomenta() const { return *m_momenta; }

  /**
   * Returns the value of \f$\text{Tr}(D^2)\f$
   */
  double traceOfDiracSquared() const;

  /**
   * Returns the value of \f$\text{Tr}(D^4)\f$
   * @return
   */
  double traceOfDirac4() const;

  /**
   * getDiracMatrix returns a fully constructed matrix representation of the
   * Dirac operator.
   */
  arma::cx_mat getDiracMatrix() const;

  // TODO: Document
  arma::cx_mat derDirac24(const int& k, const bool& herm, double g_2) const;
  // TODO: Document
  arma::cx_mat derDirac2(const int& k) const;
  // TODO: Document
  arma::cx_mat derDirac4(const int& k, const bool& herm) const;
  // TODO: Document
  std::vector<arma::cx_double>& getOmegaTable4() const { return *m_omega_table_4; }
  // TODO: Document
  void printOmegaTable4() const;
  // TODO: Document
  void randomiseMatrices(const IRng& rng_engine) const;

private:
  // CONSTANTS
  // The dimension of the H and L matrices.
  int m_dim;
  // number of H and L matrices and total number of matrices.
  int m_num_matrices, m_num_herm, m_num_antiherm;
  // size of gamma matrices
  int m_gamma_dim;

  // MATRICES
  // H and L matrices (all hermitian)
  std::unique_ptr<std::vector<arma::cx_mat>> m_matrices;
  // conjugate momenta
  std::unique_ptr<std::vector<arma::cx_mat>> m_momenta;
  // omega matrices (all hermitian)
  std::unique_ptr<std::vector<arma::cx_mat>> m_omegas;
  // epsilon: +1 for H, -1 for L
  std::unique_ptr<std::vector<int>> m_epsilons;
  // omega 4-product table
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
