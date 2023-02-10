//
// Created by Paul Druce on 10/02/2023.
//

#ifndef RFL_GSLRNG_HPP
#define RFL_GSLRNG_HPP

#include "IRng.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * This class is a thin wrap around the GSL random number generators to be used
 * in this project.
 *
 * The reason for this, is in case we want to use a different random number generator in the future
 * or if other users of this library want to switch this out for another rng.
 *
 * To do so, a new class that implements the interface defined by IRng needs to be constructed.
 */
class GslRng : public IRng {
public:
  GslRng() {
    m_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
  }
  ~GslRng() {
    if (m_rng) {
      gsl_rng_free(m_rng);
    }
  }

  double getGaussian(double sigma) const override {
    return gsl_ran_gaussian(m_rng, 1.);
  }
  double getUniform() const override {
    return gsl_rng_uniform(m_rng);
  }

private:
  gsl_rng* m_rng;
};

#endif //RFL_GSLRNG_HPP
