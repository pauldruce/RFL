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
 * @class GslRng
 *
 * @brief GSL implementation of the IRng interface.
 *
 * Wraps the GNU Scientific Library (GSL) random number generator engine.
 */
class GslRng final : public IRng {
public:
  /**
   * Constructs a GSL random number generator seeded with the current time.
   */
  GslRng() : m_rng(nullptr) {
    gsl_rng_env_setup();
    m_rng = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(m_rng, time(nullptr));
  }

  /**
   * Constructs a GSL random number generator seeded with the specified seed.
   *
   * @param seed Unsigned integer seed value.
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

  /**
   * Generates a Gaussian random variable with mean zero and standard deviation sigma.
   *
   * @param sigma Standard deviation of the Gaussian distribution.
   * @return Gaussian random value.
   */
  double getGaussian(const double sigma) const override {
    return gsl_ran_gaussian(m_rng, sigma);
  }

  /**
   * Generates a uniform random variable in the range [0, 1).
   *
   * @return Uniform random value in [0, 1).
   */
  double getUniform() const override {
    return gsl_rng_uniform(m_rng);
  }

private:
  gsl_rng* m_rng;
};

#endif//RFL_GSLRNG_HPP
