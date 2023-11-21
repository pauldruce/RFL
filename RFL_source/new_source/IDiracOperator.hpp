//
// Created by Paul Druce on 16/11/2023.
//

#ifndef IDIRACOPERATOR_HPP
#define IDIRACOPERATOR_HPP

#include "IRng.hpp"
#include <armadillo>

class IDiracOperator {
public:
  virtual ~IDiracOperator() = default;

  // TODO: Document
  virtual void randomiseMatrices(const IRng& rng_engine) const = 0;

  /**
     * Returns the Clifford type of the Dirac operator - encoded as a pair of integers
     * representing the (p,q) values of the underlying Clifford module.
     *
     * @return std::pair<int,int> of values (p,q)
     */
  virtual std::pair<int, int> getType() const = 0;

  /**
     * getMatrixDimension returns the dimension of the H and L matrices of the
     * Dirac operator
     */
  virtual int getMatrixDimension() const = 0;

  /**
     * getGammaDimension returns the dimension of the gamma matrices of the Dirac
     * operator
     */
  virtual int getGammaDimension() const = 0;

  /**
     * getNumMatrices returns the total number of H and L matrices
     */
  virtual int getNumMatrices() const = 0;

  /**
     * getNumHermitianMatrices returns the number of H matrices
     */
  virtual int getNumHermitianMatrices() const = 0;

  /**
     * getNumHermitianMatrices returns the number of L matrices
     */
  virtual int getNumAntiHermitianMatrices() const = 0;

  /**
     * This method returns a reference to the vector of H and L matrices of the Dirac operator
     */
  virtual std::vector<arma::cx_mat>& getMatrices() const = 0;

  /**
     * getEpsilons returns a reference to a vector of +/-1. The sign of the entry relates
     * to the matrix with the same index returned from getMatrices. If the value is +1 then
     * the associated matrix is Hermitian. If the value is -1, then the associated matrix
     * is anti-Hermitian.
     */
  virtual std::vector<int>& getEpsilons() const = 0;

  /**
     * getMomenta returns a reference to a vector of matrices which correspond to..
     */
  // TODO: This should probably be moved to the Hamiltonian class, or needs good documentation about it's use
  virtual std::vector<arma::cx_mat>& getMomenta() const = 0;

  /**
     * Returns the value of \f$\text{Tr}(D^2)\f$
     */
  virtual double traceOfDiracSquared() const = 0;

  /**
     * Returns the value of \f$\text{Tr}(D^4)\f$
     * @return
     */
  virtual double traceOfDirac4() const = 0;

  /**
     * getDiracMatrix returns a fully constructed matrix representation of the
     * Dirac operator.
     */
  virtual arma::cx_mat getDiracMatrix() const = 0;

  /**
     * getEigenvalues returns an armadillo vector of the eigenvalues of the Dirac Operator
     * when expressed in it's fully assembled matrix form.
     */
  virtual arma::vec getEigenvalues() const = 0;

  /**
     * Returns a vector of the Hermitian matrices, H_i, in the decomposition of
     * the Dirac Operator.
     */
  virtual std::vector<arma::cx_mat> getHermitianMatrices() const = 0;

  /**
     * Returns a vector of the Hermitian matrices, H_i, in the decomposition of
     * the Dirac Operator.
     */
  virtual std::vector<arma::cx_mat> getAntiHermitianMatrices() const = 0;

  // TODO: Document what this should do.
  // This feels like an implementation detail, but some other methods require an access to this
  // for the moment.
  virtual std::vector<arma::cx_double>& getOmegaTable4() const = 0;
};
#endif//IDIRACOPERATOR_HPP
