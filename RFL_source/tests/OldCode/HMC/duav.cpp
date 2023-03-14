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
  gsl_rng_set(engine, time(NULL));

  // create geometry
  cout << "insert p, q, dim, g2" << endl;
  Geom24 G(cin);
  G.shuffle(engine);

  double tgt = 0.6;
  double dt = 0.000001;
  cout << "initial dt: " << dt << endl;
  G.HMC_duav(100, dt, 1000, engine, tgt);
  cout << "dual averaged dt: " << dt << endl;
  double ar = G.HMC(100, dt, 100, engine);
  cout << "target acceptance rate: " << tgt << endl;
  cout << "acceptance rate: " << ar << endl;

  return 0;
}
