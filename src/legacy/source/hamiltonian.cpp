#include "Geom24.hpp"
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;
using namespace arma;

// PD: Moved into Hamiltonian class

void Geom24::sample_mom(gsl_rng* engine) {

  for (int i = 0; i < nHL; ++i) {
    // loop on indices
    for (int j = 0; j < dim; ++j) {
      double x;
      x = gsl_ran_gaussian(engine, 1.);
      mom[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < dim; ++k) {
        double a, b;
        a = gsl_ran_gaussian(engine, 1.);
        b = gsl_ran_gaussian(engine, 1.);
        mom[i](j, k) = cx_double(a, b) / sqrt(2.);
        mom[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }
}

double Geom24::calculate_K() const {
  double res = 0;

  for (int i = 0; i < nHL; ++i)
    res += trace(mom[i] * mom[i]).real();

  return res / 2;
}

double Geom24::calculate_H() const {
  return calculate_S() + calculate_K();
}
