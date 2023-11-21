//
// Created by Paul Druce on 12/05/2023.
//

#include "DiracOperator.hpp"
#include "IAlgorithm.hpp"
#include "Simulation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

class MockAlgorithm final : public IAlgorithm {
public:
  MOCK_METHOD(double, updateDirac, (const IDiracOperator& dirac), (const override));
};

TEST(SimulationTests, ConstructorDoesNotThrow) {
  ASSERT_NO_THROW(
      Simulation simulation(
          std::make_unique<DiracOperator>(1, 1, 5),
          std::make_unique<MockAlgorithm>()););
}

TEST(SimulationTests, RunCallsUpdateDirac) {
  auto algo = std::make_unique<MockAlgorithm>();
  EXPECT_CALL(*algo, updateDirac(testing::_))
      .Times(1);

  const Simulation simulation(
      std::make_unique<DiracOperator>(1, 1, 5),
      std::move(algo));

  simulation.run();
}

TEST(SimulationTests, GetDiracReturnsSameDirac) {
  auto dirac = DiracOperator(1, 1, 5);
  auto dirac_ptr = std::make_unique<DiracOperator>(dirac);
  const Simulation simulation(
      std::move(dirac_ptr),
      std::make_unique<MockAlgorithm>());

  const DiracOperator& retrieved_dirac = simulation.getDiracOperator();
  ASSERT_EQ(retrieved_dirac.getType(), std::pair(1, 1));
  ASSERT_TRUE(arma::approx_equal(
      retrieved_dirac.getDiracMatrix(),
      dirac.getDiracMatrix(),
      "absdiff",
      1e-10));
}
