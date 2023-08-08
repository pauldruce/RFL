//
// Created by Paul Druce on 12/05/2023.
//

#include "DiracOperator.hpp"
#include "IAlgorithm.hpp"
#include "Simulation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

class MockAlgorithm : public IAlgorithm {
public:
  MOCK_METHOD(double, updateDirac, (const DiracOperator& dirac), (const override));
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

  Simulation simulation(
      std::make_unique<DiracOperator>(1, 1, 5),
      std::move(algo));

  simulation.run();
}

TEST(SimulationTests, GetDiracReturnsSameDirac) {
  Simulation simulation(
      std::make_unique<DiracOperator>(1, 1, 5),
      std::make_unique<MockAlgorithm>());

  const auto& dirac = simulation.getDiracOperator();
  ASSERT_EQ(dirac.getType(), std::pair(1, 1));
}
