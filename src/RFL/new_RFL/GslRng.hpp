//
// Created by Paul Druce on 10/02/2023.
//

#ifndef RFL_GSLRNG_HPP
#define RFL_GSLRNG_HPP

#include "IRng.hpp"
#include <ctime>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/**
 * This class is a thin wrap around the GSL random number generators to be used
 * in this project.
 *
 * The reason for this, is in case we want to use a different random number generator in the future
 * or if other users of this library want to switch this out for another rng.
 *
 * To do so, a new class that implements the interface defined by IRng needs to be constructed.
 */
class GslRng final : public IRng {
public:
  /**
   * This default constructor creates a GSL random number generator that is
   * seeded by the current time.
   */
  GslRng() : m_rng(nullptr) {
    gsl_rng_env_setup();
    m_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(m_rng, time(nullptr));
  }

  /**
   * This constructor creates a GSL random number generator that is seeded by
   * the provided parameter 'seed'.
   */
  explicit GslRng(const unsigned long seed) : m_rng(nullptr) {
    m_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(m_rng, seed);
  }

  ~GslRng() override {
    if (m_rng) {
      gsl_rng_free(m_rng);
    }
  }

  double getGaussian(const double sigma) const override {
    return gsl_ran_gaussian(m_rng, sigma);
  }
  double getUniform() const override {
    return gsl_rng_uniform(m_rng);
  }

private:
  gsl_rng* m_rng;
};

#endif//RFL_GSLRNG_HPP
