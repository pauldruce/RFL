#include "Geom24.hpp"

// PD: Moved into it's own class for now.

using namespace std;
using namespace arma;

//  Standard leapfrog:
//  K/2
//  for (Nt-1) times:
//    S
//    K
//  S
//  K/2
void Geom24::leapfrog(const int& Nt, const double& dt) {
  for (int i = 0; i < nHL; ++i) {
    mat[i] += (dt / 2.) * mom[i];

    for (int j = 0; j < Nt - 1; ++j) {
      mom[i] += -dt * der_dirac24(i, true);
      mat[i] += dt * mom[i];
    }

    mom[i] += -dt * der_dirac24(i, true);
    mat[i] += (dt / 2.) * mom[i];
  }
}

void Geom24::omelyan(const int& Nt, const double& dt) {
  double xi = 0.1931833;

  for (int i = 0; i < nHL; ++i) {
    mat[i] += xi * dt * mom[i];

    for (int j = 0; j < Nt - 1; ++j) {
      mom[i] += -(dt / 2.) * der_dirac24(i, true);
      mat[i] += (1 - 2 * xi) * dt * mom[i];
      mom[i] += -(dt / 2.) * der_dirac24(i, true);
      mat[i] += 2 * xi * dt * mom[i];
    }

    mom[i] += -(dt / 2.) * der_dirac24(i, true);
    mat[i] += (1 - 2 * xi) * dt * mom[i];
    mom[i] += -(dt / 2.) * der_dirac24(i, true);
    mat[i] += xi * dt * mom[i];
  }
}
