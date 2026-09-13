#include "Cliff.hpp"
#include "Geom24.hpp"
#include <cmath>
#include <gsl/gsl_rng.h>

using namespace std;
using namespace arma;
// PD: All functions moved into Metropolis
double Geom24::MMC_duav_core(const double& scale, gsl_rng* engine, const double* s_i, double* s_f) {
  // acceptance probability
  double e;

  // metropolis
  int x = nHL * gsl_rng_uniform(engine);
  int I = dim * gsl_rng_uniform(engine);
  int J = dim * gsl_rng_uniform(engine);
  double re;
  double im = 0;
  cx_double z;
  if (I != J) {
    re = scale * (-1. + 2. * gsl_rng_uniform(engine));
    im = scale * (-1. + 2. * gsl_rng_uniform(engine));
    z = cx_double(re, im);
  } else {
    re = scale * (-1. + 2. * gsl_rng_uniform(engine));
    z = cx_double(re, 0);
  }

  double dS2 = delta2(x, I, J, z);
  double dS4 = delta4(x, I, J, z);
  double dS = g2 * dS2 + dS4;

  // metropolis test
  if (dS < 0) {
    // update matrix element
    if (I != J) {
      mat[x](I, J) += z;
      mat[x](J, I) += conj(z);
    } else {
      mat[x](I, I) += 2. * z;
    }

    // update action
    s_f[0] = s_i[0] + dS2;
    s_f[1] = s_i[1] + dS4;

    // move accepted
    e = 1;
  } else {
    e = exp(-dS);
    double uniform_p = gsl_rng_uniform(engine);

    if (e > uniform_p) {
      // update matrix element
      if (I != J) {
        mat[x](I, J) += z;
        mat[x](J, I) += conj(z);
      } else {
        mat[x](I, I) += 2. * z;
      }

      // update action
      s_f[0] = s_i[0] + dS2;
      s_f[1] = s_i[1] + dS4;
    } else {
      s_f[0] = s_i[0];
      s_f[1] = s_i[1];
    }
  }

  return e;
}

double Geom24::MMC_core(const double& scale, gsl_rng* engine, const double* s_i, double* s_f) {
  // acceptance probability
  double ret = 0;

  // metropolis
  int x = nHL * gsl_rng_uniform(engine);
  int I = dim * gsl_rng_uniform(engine);
  int J = dim * gsl_rng_uniform(engine);
  double re = 0;
  double im = 0;
  cx_double z;
  if (I != J) {
    re = scale * (-1. + 2. * gsl_rng_uniform(engine));
    im = scale * (-1. + 2. * gsl_rng_uniform(engine));
    z = cx_double(re, im);
  } else {
    re = scale * (-1. + 2. * gsl_rng_uniform(engine));
    z = cx_double(re, 0);
  }

  double dS2 = delta2(x, I, J, z);
  double dS4 = delta4(x, I, J, z);
  double dS = g2 * dS2 + dS4;

  // metropolis test
  if (dS < 0) {
    // update matrix element
    if (I != J) {
      mat[x](I, J) += z;
      mat[x](J, I) += conj(z);
    } else {
      mat[x](I, I) += 2. * z;
    }

    // update action
    s_f[0] = s_i[0] + dS2;
    s_f[1] = s_i[1] + dS4;

    // move accepted
    ret = 1;
  } else {
    double e = exp(-dS);
    double uniform_p = gsl_rng_uniform(engine);

    if (e > uniform_p) {
      // update matrix element
      if (I != J) {
        mat[x](I, J) += z;
        mat[x](J, I) += conj(z);
      } else {
        mat[x](I, I) += 2. * z;
      }

      // update action
      s_f[0] = s_i[0] + dS2;
      s_f[1] = s_i[1] + dS4;

      // move accepted
      ret = 1;
    } else {
      s_f[0] = s_i[0];
      s_f[1] = s_i[1];
    }
  }

  return ret;
}
