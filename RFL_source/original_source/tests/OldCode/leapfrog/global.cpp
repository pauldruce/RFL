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

#define N 10000
#define M 5

int main() {
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(NULL));

  // create geometry
  cout << "insert p, q, dim, g2" << endl;
  Geom24 G(cin);

  ofstream out, out1;
  string name = "GLOB";
  out.open(name + ".txt");
  out1.open(name + ".log");
  out.precision(16);
  out1.precision(16);

  double log_Nt_min = 0.7;
  double log_Nt_max = 6.2;
  vec Nt = linspace(log_Nt_min, log_Nt_max, 10);

  for (int i = 0; i < 10; ++i) {
    Nt(i) = int(exp(Nt(i)));

    double L = 0.01;
    double dt = L / Nt(i);
    cout << "Nt: " << Nt(i) << endl;
    cout << "dt: " << dt << endl;

    vec O(100);
    for (int j = 0; j < 100; ++j) {
      G.shuffle(engine);

      double Hi = G.calculate_H();
      G.omelyan(Nt(i), dt);
      double Hf = G.calculate_H();

      O(j) = fabs(Hf - Hi);
      out1 << O(j) << endl;
    }
    cout << "avg |dH|: " << mean(O) << endl;
    out << log(Nt(i)) << " " << log(mean(O)) << endl;
  }

  out.close();
  out1.close();

  return 0;
}
