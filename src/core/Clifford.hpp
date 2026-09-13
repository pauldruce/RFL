#ifndef RFL_CLIFFORD_HPP
#define RFL_CLIFFORD_HPP

#include <armadillo>
#include <vector>

/**
 * @class Clifford
 *
 * @brief Represents a Clifford module for a given signature (p, q).
 *
 * Holds gamma matrices and the chirality operator for the Clifford algebra Cl(p,q).
 */
class Clifford {
public:
  /**
   * Constructs a Clifford module of signature (p, q).
   *
   * @param p Number of Hermitian gamma matrices.
   * @param q Number of anti-Hermitian gamma matrices.
   */
  Clifford(int p, int q);

  /**
   * Returns p from the signature (p, q).
   *
   * @return Number of Hermitian gamma matrices.
   */
  int getP() const { return m_p; }

  /**
   * Returns q from the signature (p, q).
   *
   * @return Number of anti-Hermitian gamma matrices.
   */
  int getQ() const { return m_q; }

  /**
   * Returns the dimension of the generated gamma matrices.
   *
   * @return Gamma matrix dimension.
   */
  int getGammaDimension() const { return m_dim_gamma; }

  /**
   * Returns the vector of gamma matrices.
   *
   * @return Vector of gamma matrices.
   */
  std::vector<arma::cx_mat> getGammaMatrices() const { return m_gammas; }

  /**
   * Returns the gamma matrix at index i.
   *
   * @param i Zero-based index of the gamma matrix.
   * @return Gamma matrix at the specified index.
   */
  arma::cx_mat getGammaAtIndex(const int i) const { return m_gammas.at(i); }

  /**
   * Returns the chirality operator generated from the gamma matrices.
   *
   * @return Chirality operator matrix.
   */
  arma::cx_mat getChiral() const { return m_chiral; }

  /**
   * Sorts the gamma matrices by Hermiticity.
   *
   * Places Hermitian matrices first, followed by anti-Hermitian matrices.
   */
  void sortGammas();

  /**
   * Sets the maximum supported mode (p + q) for Clifford generation.
   *
   * @param max_mode Maximum mode value.
   */
  static void setMaxMode(int max_mode);

  /**
   * Returns the maximum supported mode (p + q) for Clifford generation.
   *
   * @return Maximum mode value.
   */
  static int getMaxMode();

  /**
   * Constructs a Clifford module for a base mode identifier.
   *
   * @param mode Base mode identifier (0 to 4).
   */
  explicit Clifford(int mode);

  /**
   * Copy constructor for Clifford.
   *
   * @param clifford_to_copy Clifford module to copy.
   */
  Clifford(const Clifford& clifford_to_copy);

  /**
   * Copy assignment operator for Clifford.
   *
   * @param clifford_to_copy Clifford module to copy.
   * @return Reference to this object.
   */
  Clifford& operator=(const Clifford& clifford_to_copy);

  /**
   * Multiplies this Clifford module by another using tensor products.
   *
   * @param clifford_2 Clifford module to multiply with.
   * @return Reference to this object.
   */
  Clifford& operator*=(const Clifford& clifford_2);

  /**
   * Computes the tensor product of two Clifford modules.
   *
   * @param c_1 First Clifford module.
   * @param c_2 Second Clifford module.
   * @return Combined Clifford module.
   */
  friend Clifford operator*(Clifford c_1, const Clifford& c_2) {
    c_1 *= c_2;
    return c_1;
  }
  ~Clifford() = default;

private:
  int m_p;
  int m_q;
  int m_dim_gamma;
  std::vector<arma::cx_mat> m_gammas;
  arma::cx_mat m_chiral;
  void initGammas();
};

std::ostream& operator<<(std::ostream& out, const Clifford& clifford);

/**
 * Compares two matrices to sort Hermitian matrices before anti-Hermitian matrices.
 *
 * @param m_1 First matrix.
 * @param m_2 Second matrix.
 * @return True if m_1 is Hermitian and m_2 is not.
 */
bool areHermitian(const arma::cx_mat& m_1, const arma::cx_mat& m_2);

#endif
