//
// Created by Paul Druce on 10/02/2023.
//

#ifndef RFL_IRNG_HPP
#define RFL_IRNG_HPP

/**
 * IRng is an abstract class that provides the interface that a random number generator must satisfy to be
 * used in this project.
 */
class IRng {
public:
  /**
  * getGaussian is a pure virtual method that needs to be implemented by a derived class.
  *
  * Any implementation should return a random double precision floating point number selected from a Gaussian distribution
  * with mean zero and standard deviation sigma
  *
  * @param sigma The standard deviation of the Gaussian distribution to select a number from.
  */
  virtual double getGaussian(double sigma) const = 0;

  /**
   * getUniform is a pure virtual method that needs to be implemented by a dervied class.
   *
   * Any implementation should return a double precision floating point number uniformly distributed in the range [0,1).
   * The range includes 0.0 but excludes 1.0.
   */
  virtual double getUniform() const = 0;

  virtual ~IRng() = default;
};

#endif//RFL_IRNG_HPP
