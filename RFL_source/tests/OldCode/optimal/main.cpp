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

int main() {
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));

  // create geometry
  cout << "insert p, q, dim, g2" << endl;
  Geom24 G(cin);
  G.shuffle(engine);

  // thermalize first
  double tgt = 0.8;
  double dt = 0.000001;
  G.HMC_duav(5, dt, 10000, engine, tgt);
  cout << "dual averaging complete" << endl;
  cout << "dual averaged dt: " << dt << endl;
  double ar = G.HMC(5, dt, 10000, engine);
  cout << "thermalization complete" << endl;
  cout << "acceptance rate: " << ar << endl;

  // set length
  const double L = 10;

  cout << "Nt vs ar/Nt:" << endl;
  for (int i = 1; i < 20; ++i) {
    int Nt = i;
    Geom24 G1 = G;
    ar = G1.HMC(Nt, L / Nt, 1000, engine);
    cout << Nt << " " << ar / Nt << endl;
  }

  return 0;
}
