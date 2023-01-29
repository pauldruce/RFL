#ifndef CLIFFORD_HPP
#define CLIFFORD_HPP

#include <armadillo>
#include <vector>

class Clifford {
 public:
  // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR
  explicit Clifford(int mode);
  Clifford(int p, int q);
  Clifford(const Clifford &C);
  Clifford &operator=(const Clifford &C);
  ~Clifford() = default;
  // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR

  // ============== OPERATORS
  Clifford &operator*=(const Clifford &C);
  friend Clifford operator*(Clifford C1, const Clifford &C2) {
	C1 *= C2;
	return C1;
  }
  // ============== OPERATORS

  // ============== GET METHODS
  int get_p() const { return p; }
  int get_q() const { return q; }
  int get_dim_gamma() const { return dim_gamma; }
  std::vector<arma::cx_mat> get_gammas() const { return gammas; }
  arma::cx_mat get_gamma(int i) const { return gammas.at(i); }
  arma::cx_mat get_chiral() const { return chiral; }
  // ============== GET METHODS


  // ============== OTHER METHODS
  void sort_gamma();
  // ============== OTHER METHODS

 private:
  int p;
  int q;
  int dim_gamma;
  std::vector<arma::cx_mat> gammas;
  arma::cx_mat chiral;
  void init_gamma();

};

std::ostream &operator<<(std::ostream &out, const Clifford &C);
//void decomp(int p, int q, int &dec);
bool hermiticity(const arma::cx_mat &M1, const arma::cx_mat &M2);

#endif

