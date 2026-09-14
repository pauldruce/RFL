#ifndef RFL_CLIFFORD_HPP
#define RFL_CLIFFORD_HPP

#include <armadillo>
#include <vector>

class Cliff {
public:
  // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR
  explicit Cliff(int mode);

  Cliff(int p, int q);

  Cliff(const Cliff& C);

  Cliff& operator=(const Cliff& C);

  ~Cliff() = default;
  // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR

  // ============== OPERATORS
  Cliff& operator*=(const Cliff& C);

  friend Cliff operator*(Cliff C1, const Cliff& C2) {
    C1 *= C2;
    return C1;
  }
  // ============== OPERATORS

  // ============== GET METHODS
  int get_p() const { return p; }

  int get_q() const { return q; }

  int get_dim_gamma() const { return dim_gamma; }

  arma::cx_mat get_gamma(int i) const { return gamma.at(i); }

  std::vector<arma::cx_mat> get_gamma() const { return gamma; }

  arma::cx_mat get_chiral() const { return chiral; }
  // ============== GET METHODS

  // ============== OTHER METHODS
  void sort_gamma();
  // ============== OTHER METHODS

private:
  int p;
  int q;

  int dim_gamma;

  std::vector<arma::cx_mat> gamma;
  arma::cx_mat chiral;

  void init_gamma();
};

std::ostream& operator<<(std::ostream& out, const Cliff& C);

//void decomp(int p, int q, int &dec);

//static bool areHermitian(const arma::cx_mat& m_1, const arma::cx_mat& m_2);

#endif
