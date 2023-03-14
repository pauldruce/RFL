//
// Created by Paul Druce on 26/06/2022.
//

#include "Geom24.hpp"
#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

static auto engine = gsl_rng_alloc(gsl_rng_ranlxd1);

static void CompareActions(int p, int q, int dim, double g2) {
  constexpr int N = 100;
  Geom24 G(p, q, dim, g2);
  for (int i = 0; i < N; ++i) {
    G.shuffle(engine);
    double d = G.get_dim();

    double S1 = G.calculate_S() / (d * d);
    double S2 = G.calculate_S_from_dirac() / (d * d);
    EXPECT_TRUE(fabs(S1 - S2) < 1e-8) << "Methods differ more then 1e-8";
  }
}

TEST(ActionTests, DifferenceIsSmall) {
  typedef struct ActionParameters {
    int p;
    int q;
    int dim;
    double g2;
  } ActionParameters;

  const ActionParameters testing_params[] =
      {
          {1, 1, 5, -3.0},
          {2, 2, 6, -2.2},
          {1, 2, 7, -0.5},
          {2, 1, 4, -0.1},
          {0, 5, 4, -2.8}};

  for (auto& d : testing_params) {
    CompareActions(d.p, d.q, d.dim, d.g2);
  }
}

//! PARAMETERISED TESTING USING GTEST
//
//typedef struct ActionParameters
//{
//   int p;
//   int q;
//   int dim;
//   double g2;
//} ActionParameters;
//
//class ActionTests : public ::testing::TestWithParam<ActionParameters>
//{
//};
//
//TEST_P(ActionTests, DifferenceIsSmall)
//{
//   ActionParameters action_params = GetParam();
//   CompareActions(action_params.p, action_params.q, action_params.dim, action_params.g2);
//}
//
//const ActionParameters testing_params[] =
//    {
//        {1, 1, 5, -3.0},
//        {2, 2, 6, -2.2},
//        {1, 2, 7, -0.5},
//        {2, 1, 4, -0.1},
//        {0, 5, 4, -2.8}};
//
//INSTANTIATE_TEST_SUITE_P(ActionTestInstantiated,
//                         ActionTests,
//                         testing::ValuesIn(testing_params));