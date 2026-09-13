#include "Geom24.hpp"
#include <armadillo>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <iostream>

using namespace std;
using namespace arma;

#define N 100

int main() {
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));

  // create geometry
  cout << "insert p, q, dim, g2" << endl;
  Geom24 G(cin);
  cout.precision(16);

  G.shuffle(engine);
  Geom24 G_old = G;

  double Si = G.calculate_S();
  double Ki = G.calculate_K();
  G.leapfrog(100, 0.00001);
  G.reverse_mom();
  G.leapfrog(100, 0.00001);
  G.reverse_mom();
  double Sf = G.calculate_S();
  double Kf = G.calculate_K();
  cout.precision(16);
  cout << "initial potential|kinetic|energy: " << Si << "|" << Ki << "|" << Si + Ki << endl;
  cout << "final   potential|kinetic|energy: " << Sf << "|" << Kf << "|" << Sf + Kf << endl;
  for (int i = 0; i < G.get_nHL(); ++i) {
    cout << "mat[" + to_string(i) + "]: " << approx_equal(G_old.get_mat(i), G.get_mat(i), "absdiff", 0.001) << endl;
    cout << "mom[" + to_string(i) + "]: " << approx_equal(G_old.get_mom(i), G.get_mom(i), "absdiff", 0.001) << endl;
  }

  return 0;
}
