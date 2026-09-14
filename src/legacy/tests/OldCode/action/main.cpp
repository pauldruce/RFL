#include "Geom24.hpp"
#include <armadillo>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <iostream>

using namespace std;
using namespace arma;

#define N 1000

int main() {
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));

  // create geometry
  cout << "insert p, q, dim, g2" << endl;
  Geom24 G(cin);

  for (int i = 0; i < N; ++i) {
    if (!(i % 100)) {
      cout << "iteration " << i << endl;
    }

    G.shuffle(engine);
    double d = G.get_dim();

    double S1 = G.calculate_S() / (d * d);
    double S2 = G.calculate_S_from_dirac() / (d * d);

    cout.precision(16);
    if (fabs(S1 - S2) > 1e-8) {
      cout << S1 << " " << S2 << " " << S1 - S2 << endl;
    }
  }

  return 0;
}
