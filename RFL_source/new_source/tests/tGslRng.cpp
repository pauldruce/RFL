//
// Created by Paul Druce on 11/02/2023.
//

#include "GslRng.hpp"
#include <gtest/gtest.h>

TEST(GslRngTests, SettingSeedProducesSameValues) {
  GslRng rng(1);

  double expected_values[5] = {0.83451879245814453,
                               0.61670202724383927,
                               0.44438336146091828,
                               0.038517618778346474,
                               0.5896974345675261};
  std::cout.precision(17);
  for (double expected_value : expected_values) {
    EXPECT_DOUBLE_EQ(rng.getUniform(), expected_value);
  }
}

TEST(GslRngTests, NoSeedProducesDifferentValues) {
  GslRng rng;

  // Get first 5 numbers
  double first_values[10];
  for (double& first_value : first_values) {
    first_value = rng.getUniform();
  }

  GslRng new_rng;
  for (double orig_val : first_values) {
    EXPECT_LE(std::abs(orig_val - new_rng.getUniform()), 1e-8);
  }
}