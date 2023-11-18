#ifndef RFL_CLIFFORD_HPP
#define RFL_CLIFFORD_HPP

#include <armadillo>
#include <vector>

/**
 * A class that represents a Clifford module of a specific type.
 */
class Clifford {
public:
  /**
   * The constructor for this class. This creates a Clifford modules of type
   * \f$(p,q)\f$.
   */
  Clifford(int p, int q);

  /**
   * getP returns the value of \f$p\f$ from the Clifford type \f$(p,q)\f$
   */
  int getP() const { return m_p; }

  /**
   * getQ returns the value of \f$q\f$ from the Clifford type \f$(p,q)\f$
   */
  int getQ() const { return m_q; }

  /**
   * getGammaDimension returns the dimension of the generated gamma matrices.
   */
  int getGammaDimension() const { return m_dim_gamma; }

  /**
   * getGammaMatrices returns a vector of matrices which is a copy of the gamma
   * matrices.
   */
  std::vector<arma::cx_mat> getGammaMatrices() const { return m_gammas; }

  /**
   * getGammaAtIndex returns a copy of a gamma matrix located in the array at index
   * i.
   */
  arma::cx_mat getGammaAtIndex(const int i) const { return m_gammas.at(i); }

  /**
   * getChiral returns a copy of the chirality matrix generated from the gamma
   * matrices.
   * See the reference "Spin Geometry" by Lawson and Michelson for more techincal
   * details.
   */
  arma::cx_mat getChiral() const { return m_chiral; }

  /**
   * sortGammas orders the gamma matrices in the vector that houses them. The
   * ordering is based on Hermiticity. The Hermitian matrices will be first,
   * with the anti-Hermitian matrices afterwards in the vector.
   */
  void sortGammas();

  explicit Clifford(int mode);
  Clifford(const Clifford& clifford_to_copy);
  Clifford& operator=(const Clifford& clifford_to_copy);
  Clifford& operator*=(const Clifford& clifford_2);
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
bool areHermitian(const arma::cx_mat& m_1, const arma::cx_mat& m_2);

#endif
