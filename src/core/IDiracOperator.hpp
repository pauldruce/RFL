//
// Created by Paul Druce on 16/11/2023.
//

#ifndef IDIRACOPERATOR_HPP
#define IDIRACOPERATOR_HPP

#include "IRng.hpp"
#include <armadillo>

/**
 * @interface IDiracOperator
 *
 * @brief Interface for finite Dirac operators.
 *
 * Defines the contract for Dirac operators in finite noncommutative geometries.
 */
class IDiracOperator {
public:
  virtual ~IDiracOperator() = default;

  /**
   * Randomises all internal matrix elements using the provided random number generator.
   *
   * @param rng_engine Random number generator engine.
   */
  virtual void randomiseMatrices(const IRng& rng_engine) const = 0;

  /**
   * Returns the signature of the Dirac operator as a pair (p,q).
   *
   * @return std::pair<int,int> containing (p,q).
   */
  virtual std::pair<int, int> getType() const = 0;

  /**
   * Returns the matrix dimension of the H and L matrices.
   *
   * @return Matrix dimension.
   */
  virtual int getMatrixDimension() const = 0;

  /**
   * Returns the dimension of the gamma matrices.
   *
   * @return Gamma matrix dimension.
   */
  virtual int getGammaDimension() const = 0;

  /**
   * Returns the total number of H and L matrices.
   *
   * @return Number of matrices.
   */
  virtual int getNumMatrices() const = 0;

  /**
   * Returns the number of Hermitian matrices (H).
   *
   * @return Number of Hermitian matrices.
   */
  virtual int getNumHermitianMatrices() const = 0;

  /**
   * Returns the number of anti-Hermitian matrices (L).
   *
   * @return Number of anti-Hermitian matrices.
   */
  virtual int getNumAntiHermitianMatrices() const = 0;

  /**
   * Returns a reference to the vector of H and L matrices.
   *
   * @return Reference to vector of matrices.
   */
  virtual std::vector<arma::cx_mat>& getMatrices() const = 0;

  /**
   * Returns a reference to the vector of signs (+1 or -1).
   *
   * The sign matches the matrix at the same index in getMatrices().
   * A value of +1 indicates a Hermitian matrix (H).
   * A value of -1 indicates an anti-Hermitian matrix (L).
   *
   * @return Reference to vector of signs.
   */
  virtual std::vector<int>& getEpsilons() const = 0;

  /**
   * Returns a reference to the vector of conjugate momenta matrices.
   *
   * @return Reference to vector of momenta matrices.
   */
  // TODO: Move to Hamiltonian class or document usage.
  virtual std::vector<arma::cx_mat>& getMomenta() const = 0;

  /**
   * Calculates the trace of the Dirac operator squared, \f$\text{Tr}(D^2)\f$.
   *
   * @return Trace value as a double.
   */
  virtual double traceOfDiracSquared() const = 0;

  /**
   * Calculates the trace of the fourth power of the Dirac operator, \f$\text{Tr}(D^4)\f$.
   *
   * @return Trace value as a double.
   */
  virtual double traceOfDirac4() const = 0;

  /**
   * Constructs and returns the assembled matrix representation of the Dirac operator.
   *
   * @return Assembled Dirac operator matrix.
   */
  virtual arma::cx_mat getDiracMatrix() const = 0;

  /**
   * Calculates and returns the eigenvalues of the assembled Dirac operator matrix.
   *
   * @return Armadillo vector of real eigenvalues.
   */
  virtual arma::vec getEigenvalues() const = 0;

  /**
   * Returns a vector of the Hermitian matrices (H_i) in the Dirac operator.
   *
   * @return Vector of Hermitian matrices.
   */
  virtual std::vector<arma::cx_mat> getHermitianMatrices() const = 0;

  /**
   * Returns a vector of the anti-Hermitian matrices (L_j) in the Dirac operator.
   *
   * @return Vector of anti-Hermitian matrices.
   */
  virtual std::vector<arma::cx_mat> getAntiHermitianMatrices() const = 0;

  /**
   * Returns a reference to the four-product table of omega matrices.
   *
   * @return Reference to four-product table vector.
   */
  virtual std::vector<arma::cx_double>& getOmegaTable4() const = 0;
};
#endif//IDIRACOPERATOR_HPP
