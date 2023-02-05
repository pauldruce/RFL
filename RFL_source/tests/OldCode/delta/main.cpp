#include "Geom24.hpp"
#include <armadillo>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace arma;

#define N 1000

int main() {
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(NULL));

  // create geometry
  cout << "insert p, q, dim, g2" << endl;
  Geom24 G(cin);

  for (int i = 0; i < N; ++i) {
    if (!(i % 100)) {
      cout << "iteration " << i << endl;
    }

    G.shuffle(engine);
    double Si = G.calculate_S();

    // random entry
    int x = G.get_nHL() * gsl_rng_uniform(engine);
    int I = G.get_dim() * gsl_rng_uniform(engine);
    int J = G.get_dim() * gsl_rng_uniform(engine);
    double re = 0;
    double im = 0;
    cx_double z;
    if (I != J) {
      re = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
      im = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
      z = cx_double(re, im);
    } else {
      re = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
      z = cx_double(re, 0);
    }

    double dS = G.delta24(x, I, J, z);

    cx_mat m = G.get_mat(x);
    m(I, J) += z;
    m(J, I) += conj(z);
    G.set_mat(x, m);

    double Sf = G.calculate_S();

    cout.precision(16);
    if (fabs(Sf - Si - dS) > 1e-8) {
      cout << Sf - Si << " " << dS << " " << Sf - Si - dS << endl;
    }
  }

  return 0;
}
