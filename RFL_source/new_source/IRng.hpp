//
// Created by Paul Druce on 10/02/2023.
//

#ifndef RFL_IRNG_HPP
#define RFL_IRNG_HPP

/**
 * This abstract class provides an interface that a random number generator must satisfy to be
 * used in this project.
 */
class IRng {
public:
  /**
   * This function returns a Gaussian random variate, with mean zero and standard deviation sigma.
   */
  virtual double getGaussian(double sigma) const = 0;

  /**
   * This function returns a double precision floating point number uniformly distributed in the range [0,1).
   * The range includes 0.0 but excludes 1.0
   */
  virtual double getUniform() const = 0;

  virtual ~IRng() = default;
};

#endif //RFL_IRNG_HPP
