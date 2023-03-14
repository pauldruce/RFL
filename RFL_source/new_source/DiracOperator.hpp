//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_DIRACOPERATOR_HPP
#define RFL_DIRACOPERATOR_HPP

#include "Clifford.hpp"
#include "IRng.hpp"
#include <armadillo>

class DiracOperator {
public:
  // CONSTRUCTOR/DESTRUCTORS
  DiracOperator() = delete;
  DiracOperator(int p, int q, int dim);
  ~DiracOperator();

public:
  int getMatrixDimension() const {return m_dim;};
  int getGammaDimension() const {return m_gamma_dim;};
  int getNumMatrices() const {return m_num_matrices;} ;
  int getNumHermitianMatrices() const {return m_num_herm;};
  int getNumAntiHermitianMatrices() const {return m_num_antiherm;};

  // METHODS
  arma::cx_mat* getMatrices() const { return m_matrices; } // H and L matrices in the Dirac Operator
  int* getEpsilons() const { return m_epsilons; }           //
  arma::cx_mat* getMomenta() const { return m_momenta; } // Conjugate momenta for use with Hamiltonian method

  arma::cx_mat getDiracMatrix() const;

  arma::cx_mat derDirac24(const int& k, const bool& herm, double g_2) const;
  arma::cx_mat derDirac2(const int& k) const;
  arma::cx_mat derDirac4(const int& k, const bool& herm) const;

  arma::cx_double* getOmegaTable4() const { return m_omega_table_4; }
  void printOmegaTable4() const;

  void randomiseMatrices(const IRng& rng_engine);

private:
  // CONSTANTS
  int m_dim;         // The dimension of the H and L matrices.
  int m_num_matrices, m_num_herm, m_num_antiherm; // number of H and L matrices and total number of matrices.
  int m_gamma_dim;   // size of gamma matrices

  // MATRICES
  arma::cx_mat* m_matrices;   // H and L matrices (all hermitian)
  arma::cx_mat* m_momenta;   // conjugate momenta
  arma::cx_mat* m_omegas; // omega matrices (all hermitian)

  int* m_epsilons;                         // epsilon: +1 for H, -1 for L
  arma::cx_double* m_omega_table_4{}; // omega 4-product table

  void initOmegaTable4();

  arma::cx_mat computeB4(const int& k,
                         const int& i_2,
                         const int& i_3,
                         const int& i_4,
                         const double& cliff,
                         const bool& neg) const;

  arma::cx_mat computeB2(const int& k, const int& i) const;
  arma::cx_mat computeB(const int& k) const;
};
#endif // RFL_DIRACOPERATOR_HPP
