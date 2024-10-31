//
// Created by Paul Druce on 27/12/2022.
//
#include "BarrettGlaser/Action.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"
#include <gtest/gtest.h>

static void CompareActions(const int p, const int q, const int dim, const double g_2) {
  constexpr int num_of_test_repeats = 100;
  const Action action(g_2);
  const DiracOperator dirac(p, q, dim);

  for (int i = 0; i < num_of_test_repeats; ++i) {
    const GslRng rng;
    dirac.randomiseMatrices(rng);
    const double d_2 = dirac.getMatrixDimension() * dirac.getMatrixDimension();

    auto s_1 = action.calculateS(dirac) / d_2;

    auto s_2 = action.calculateSFromDirac(dirac) / d_2;

    EXPECT_TRUE(fabs(s_1 - s_2) < 1e-8) << "Methods differ more then 1e-8.\n"
                                        << "Optimised calculation = " << s_1 << std::endl
                                        << "Direct calculation = " << s_2 << std::endl;
  }
}

TEST(ActionTests, ActionMethodsDontDiffer) {
  typedef struct ActionParameters {
    int p;
    int q;
    int dim;
    double g_2;
  } ActionParameters;

  const ActionParameters testing_params[] =
      {
          {1, 1, 5, -3.0},
          {2, 2, 6, -2.2},
          {1, 2, 7, -0.5},
          {2, 1, 4, -0.1},
          {0, 5, 4, -2.8}};

  for (const auto& [p, q, dim, g_2] : testing_params) {
    CompareActions(p, q, dim, g_2);
  }
}

TEST(ActionTests, Set_g2) {
  constexpr double initial_g_2 = 1.1;
  Action action(initial_g_2);
  EXPECT_EQ(action.getG2(), initial_g_2);

  constexpr double new_g_2 = -3.4;
  action.setG2(new_g_2);
  EXPECT_EQ(action.getG2(), new_g_2);
}

TEST(ActionTests, Set_g4) {
  constexpr double initial_g_4 = 1.1;
  Action action(1.0, initial_g_4);
  EXPECT_EQ(action.getG4(), initial_g_4);

  constexpr double new_g_4 = -3.4;
  action.setG4(new_g_4);
  EXPECT_EQ(action.getG4(), new_g_4);
}

TEST(ActionTests, SetParameters) {
  constexpr double initial_g_2 = 1.1, initial_g_4 = 2.2;
  Action action(initial_g_2, initial_g_4);

  EXPECT_EQ(action.getG2(), initial_g_2);
  EXPECT_EQ(action.getG4(), initial_g_4);

  constexpr double new_g_2 = -5.0, new_g_4 = 7.2;
  action.setParams(new_g_2, new_g_4);

  EXPECT_EQ(action.getG2(), new_g_2);
  EXPECT_EQ(action.getG4(), new_g_4);
}

TEST(ActionTests, CreateWithNoParams) {
  EXPECT_NO_THROW(
      const Action action;
      EXPECT_EQ(action.getG2(), 0);
      EXPECT_EQ(action.getG4(), 0););
}