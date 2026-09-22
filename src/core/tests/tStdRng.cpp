//
// Created for Random Fuzzy Library (RFL).
//

#include "../StdRng.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

TEST(StdRngTests, SettingSeedProducesSameValues) {
  const StdRng rng1(42UL);
  const StdRng rng2(42UL);

  for (int i = 0; i < 100; ++i) {
    EXPECT_DOUBLE_EQ(rng1.getUniform(), rng2.getUniform());
    EXPECT_DOUBLE_EQ(rng1.getGaussian(1.5), rng2.getGaussian(1.5));
  }
}

TEST(StdRngTests, DifferentSeedsProduceDifferentValues) {
  const StdRng rng1(42UL);
  const StdRng rng2(43UL);

  bool any_different = false;
  for (int i = 0; i < 10; ++i) {
    if (rng1.getUniform() != rng2.getUniform()) {
      any_different = true;
      break;
    }
  }
  EXPECT_TRUE(any_different);
}

TEST(StdRngTests, UnseededProducesDifferentSequences) {
  const StdRng rng1;
  const StdRng rng2;

  bool any_different = false;
  for (int i = 0; i < 10; ++i) {
    if (rng1.getUniform() != rng2.getUniform()) {
      any_different = true;
      break;
    }
  }
  EXPECT_TRUE(any_different);
}

TEST(StdRngTests, UniformDistributionBounds) {
  const StdRng rng(12345UL);
  for (int i = 0; i < 1000; ++i) {
    const double val = rng.getUniform();
    EXPECT_GE(val, 0.0);
    EXPECT_LT(val, 1.0);
  }
}

TEST(StdRngTests, GaussianDistributionMoments) {
  const StdRng rng(98765UL);
  constexpr int num_samples = 10000;
  constexpr double sigma = 2.0;

  double sum = 0.0;
  double sum_sq = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    const double x = rng.getGaussian(sigma);
    sum += x;
    sum_sq += x * x;
  }

  const double mean = sum / num_samples;
  const double variance = (sum_sq / num_samples) - (mean * mean);

  // Mean should be near 0 with high statistical confidence (|mean| < 3 * sigma / sqrt(N) ≈ 0.06)
  EXPECT_NEAR(mean, 0.0, 0.1);
  // Variance should be near sigma^2 = 4.0
  EXPECT_NEAR(variance, sigma * sigma, 0.2);
}

TEST(StdRngTests, DiscreteUniformIntBounds) {
  const StdRng rng(54321UL);
  constexpr uint64_t min_val = 5;
  constexpr uint64_t max_val = 15;

  std::vector<bool> seen(max_val - min_val + 1, false);
  for (int i = 0; i < 1000; ++i) {
    const uint64_t val = rng.getUniformInt(min_val, max_val);
    EXPECT_GE(val, min_val);
    EXPECT_LE(val, max_val);
    seen[val - min_val] = true;
  }

  for (const bool hit : seen) {
    EXPECT_TRUE(hit) << "Every integer in range [min, max] should be sampled across 1000 trials.";
  }
}
