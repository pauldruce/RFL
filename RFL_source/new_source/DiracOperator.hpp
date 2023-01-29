//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_DIRACOPERATOR_HPP
#define RFL_DIRACOPERATOR_HPP

#include <armadillo>
#include "Clifford.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class DiracOperator {
 public:
  // CONSTRUCTOR/DESTRUCTORS
  DiracOperator() = delete;
  DiracOperator(int p, int q, int dim);
  ~DiracOperator();

 public:
  // CONSTANTS
  int dim; // The dimension of the H and L matrices.
  int nHL, nH, nL; // number of H and L matrices (nH and nL) and total number of matrices (nHL)
  int dim_omega; // size of gamma matrices

  // METHODS
  arma::cx_mat *get_mats() const { return mat; } // H and L matrices in the Dirac Operator
  int *get_eps() const { return eps; } //
  arma::cx_mat *get_moms() const { return mom; } // Conjugate momenta for use with Hamiltonian method

  arma::cx_mat build_dirac() const;

  arma::cx_mat der_dirac24(const int &k, const bool &herm, double g2) const;
  arma::cx_mat der_dirac2(const int &k) const;
  arma::cx_mat der_dirac4(const int &k, const bool &herm) const;

  arma::cx_double *get_omega_table_4() const { return omega_table_4; }
  void print_omega_table_4() const;

  void randomise(gsl_rng *engine);

 private:
  // MATRICES
  arma::cx_mat *mat; // H and L matrices (all hermitian)
  arma::cx_mat *mom; // conjugate momenta
  arma::cx_mat *omega; // omega matrices (all hermitian)

  int *eps; // epsilon: +1 for H, -1 for L
  arma::cx_double *omega_table_4{}; // omega 4-product table

  void init_omega_table_4();

  arma::cx_mat compute_B4(const int &k,
						  const int &i2,
						  const int &i3,
						  const int &i4,
						  const double &cliff,
						  const bool &neg) const;

  arma::cx_mat compute_B2(const int &k, const int &i) const;
  arma::cx_mat compute_B(const int &k) const;
};
#endif // RFL_DIRACOPERATOR_HPP
