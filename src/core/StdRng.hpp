//
// Created for Random Fuzzy Library (RFL).
//

#ifndef RFL_CORE_STDRNG_HPP
#define RFL_CORE_STDRNG_HPP

#include "IRng.hpp"
#include <chrono>
#include <cstdint>
#include <random>

/**
 * @class StdRng
 *
 * @brief Standard C++ implementation of the IRng interface.
 *
 * Generates uniform and Gaussian random numbers using the standard C++
 * 64-bit Mersenne Twister engine (std::mt19937_64).
 */
class StdRng final : public IRng {
public:
  /**
   * Constructs an instance seeded with non-deterministic system entropy.
   */
  StdRng() : m_engine(generateSeed()) {}

  /**
   * Constructs an instance with an explicit seed for reproducible simulations.
   *
   * @param seed 64-bit unsigned integer seed value.
   */
  explicit StdRng(const uint64_t seed) : m_engine(seed) {}

  ~StdRng() override = default;

  /**
   * Generates a Gaussian random variable with mean zero and standard deviation sigma.
   *
   * @param sigma Standard deviation of the Gaussian distribution.
   * @return Gaussian random value.
   */
  double getGaussian(const double sigma) const override {
    return sigma * m_normal(m_engine);
  }

  /**
   * Generates a uniform random variable in the half-open interval [0, 1).
   *
   * @return Uniform random value in [0, 1).
   */
  double getUniform() const override {
    return m_uniform(m_engine);
  }

  /**
   * Generates a discrete uniform integer in the closed interval [min, max].
   *
   * @param min Lower bound (inclusive).
   * @param max Upper bound (inclusive).
   * @return Uniform random integer in [min, max].
   */
  uint64_t getUniformInt(const uint64_t min, const uint64_t max) const override {
    std::uniform_int_distribution<uint64_t> dist(min, max);
    return dist(m_engine);
  }

private:
  static uint64_t generateSeed() {
    std::random_device rd;
    const uint64_t part1 = static_cast<uint64_t>(rd());
    const uint64_t part2 = static_cast<uint64_t>(rd());
    uint64_t seed = (part1 << 32) | part2;
    if (seed == 0) {
      const auto now = std::chrono::high_resolution_clock::now();
      seed = static_cast<uint64_t>(now.time_since_epoch().count());
    }
    return seed;
  }

  mutable std::mt19937_64 m_engine;
  mutable std::uniform_real_distribution<double> m_uniform{0.0, 1.0};
  mutable std::normal_distribution<double> m_normal{0.0, 1.0};
};

#endif// RFL_CORE_STDRNG_HPP
