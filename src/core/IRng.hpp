//
// Created by Paul Druce on 10/02/2023.
//

#ifndef RFL_IRNG_HPP
#define RFL_IRNG_HPP

#include <cstdint>

/**
 * @interface IRng
 *
 * @brief Interface for random number generator engines.
 *
 * Defines methods for drawing Gaussian, uniform, and discrete uniform random numbers.
 */
class IRng {
public:
  /**
   * Generates a Gaussian random variable with mean zero and standard deviation sigma.
   *
   * @param sigma Standard deviation of the Gaussian distribution.
   * @return Gaussian random value.
   */
  virtual double getGaussian(double sigma) const = 0;

  /**
   * Generates a uniform random variable in the range [0, 1).
   *
   * @return Uniform random value in [0, 1).
   */
  virtual double getUniform() const = 0;

  /**
   * Generates a discrete uniform integer in the closed interval [min, max].
   *
   * @param min Lower bound (inclusive).
   * @param max Upper bound (inclusive).
   * @return Uniform random integer in [min, max].
   */
  virtual uint64_t getUniformInt(uint64_t min, uint64_t max) const = 0;

  virtual ~IRng() = default;
};

#endif//RFL_IRNG_HPP
