#ifndef RFL_CLIFFORD_HPP
#define RFL_CLIFFORD_HPP

#include <armadillo>
#include <vector>

class Clifford {
public:
  // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR
  explicit Clifford(int mode);
  Clifford(int p, int q);
  Clifford(const Clifford& clifford_to_copy);
  Clifford& operator=(const Clifford& clifford_to_copy);
  ~Clifford() = default;
  // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR

  // ============== OPERATORS
  Clifford& operator*=(const Clifford& clifford_2);
  friend Clifford operator*(Clifford c_1, const Clifford& c_2) {
    c_1 *= c_2;
    return c_1;
  }
  // ============== OPERATORS

  // ============== GET METHODS
  int getP() const { return m_p; }
  int getQ() const { return m_q; }
  int getGammaDimension() const { return m_dim_gamma; }
  std::vector<arma::cx_mat> getGammaMatrices() const { return m_gammas; }
  arma::cx_mat getGammaAtIndex(int i) const { return m_gammas.at(i); }
  arma::cx_mat getChiral() const { return m_chiral; }
  // ============== GET METHODS

  // ============== OTHER METHODS
  void sortGammas();
  // ============== OTHER METHODS

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
