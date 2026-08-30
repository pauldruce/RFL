//
// Created by Paul Druce on 10/02/2023.
//

#ifndef RFL_IRNG_HPP
#define RFL_IRNG_HPP

/**
 * @interface IRng
 *
 * @brief Interface for random number generator engines.
 *
 * Defines methods for drawing Gaussian and uniform random numbers.
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

  virtual ~IRng() = default;
};

#endif//RFL_IRNG_HPP
