//
// Created by Paul Druce on 29/06/2022.
//

#include "Geom24.hpp"
#include "gsl/gsl_rng.h"
#include <armadillo>
#include <gtest/gtest.h>

namespace DeltaTests {
using namespace arma;

static auto engine = gsl_rng_alloc(gsl_rng_ranlxd1);
constexpr int N = 10;

static void CompareDeltaToDifference(int p, int q, int dim, double g2) {
  Geom24 G(p, q, dim, g2);

  for (int i = 0; i < N; ++i) {
    G.shuffle(engine);
    double Si = G.calculate_S();

    // random entry
    int x = (int)(G.get_nHL() * gsl_rng_uniform(engine));
    int I = (int)(G.get_dim() * gsl_rng_uniform(engine));
    int J = (int)(G.get_dim() * gsl_rng_uniform(engine));
    double re;
    double im;
    cx_double z;
    if (I != J) {
      re = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
      im = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
      z = cx_double(re, im);
    } else {
      re = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
      z = cx_double(re, 0);
    }

    double dS = G.delta24(x, I, J, z);

    cx_mat m = G.get_mat(x);
    m(I, J) += z;
    m(J, I) += conj(z);
    G.set_mat(x, m);

    double Sf = G.calculate_S();
    EXPECT_TRUE(fabs((Sf - Si) - dS) < 1e-7)
        << "     Sf-Si =" << Sf - Si << ", dS = " << dS << " and |(Sf-Si) - dS| = "
        << fabs(Sf - Si - dS);
  }
}

TEST(DeltaTest, DifferenceBetweenTwoActionMethodsIsSmall) {
  typedef struct DeltaTestParameters {
    int p;
    int q;
    int dim;
    double g2;
  } DeltaTestParameters;

  std::vector<DeltaTestParameters> test_params;

  // Create (p,q) combo's upto n = p+q=5.
  std::pair<int, int> pq_pairs[] = {
      {1, 0},
      {0, 1},
      {2, 0},
      {1, 1},
      {0, 2},
      {3, 0},
      {2, 1},
      {1, 2},
      {0, 3},
      {4, 0},
      {3, 1},
      {2, 2},
      {1, 3},
      {0, 4}};

  for (auto pair : pq_pairs) {
    for (int dim = 2; dim < 6; dim++) {
      for (double g2 = -5.0; g2 < 0.0; g2 += 1.0) {
        DeltaTestParameters params = {pair.first, pair.second, dim, g2};
        test_params.emplace_back(params);
        ASSERT_NO_THROW(CompareDeltaToDifference(pair.first, pair.second, dim, g2))
            << "An exception was thrown, most likely by constructor for Geom24";
      }
    }
  }
}
}// namespace DeltaTests

// PARAMETERISED TESTS USING GTEST
//
//   struct DeltaTestParameters
//   {
//      int p;
//      int q;
//      int dim;
//      double g2;
//   };
//
//
//   class DeltaTests : public ::testing::TestWithParam<DeltaTestParameters>
//   {
//   };
//
//   TEST_P(DeltaTests, DifferenceBetweenTwoMethodsIsSmall)
//   {
//      DeltaTestParameters params = GetParam();
//      ASSERT_NO_THROW(CompareDeltaToDifference(params.p, params.q, params.dim, params.g2))
//                           << "An exception was thrown, most likely by constructor for Geom24";
//   }
//
//   static std::vector<DeltaTestParameters> GetTestParameters()
//   {
//      std::vector<DeltaTestParameters> test_params;
//
//      // Create (p,q) combo's upto n = p+q=5.
//      std::pair<int, int> pq_pairs[] = {
//            {1, 0},{0, 1},
//            {2, 0},{1, 1},{0, 2},
//            {3, 0},{2, 1},{1, 2},{0, 3},
//            {4, 0},{3, 1},{2, 2},{1, 3},{0, 4},
//            {5, 0},{4, 1},{3, 2},{2, 3},{1, 4},{0, 5}
//      };
//
//      for (auto pair: pq_pairs)
//      {
//         for (int dim = 2; dim < 10; dim++)
//         {
//            for (double g2 = -5.0; g2 < 0.0; g2 += 1.0)
//            {
//               DeltaTestParameters params = {pair.first, pair.second, dim, g2};
//               test_params.emplace_back(params);
//            }
//         }
//      }
//      return test_params;
//   }
//
//   INSTANTIATE_TEST_SUITE_P(DeltaTestsInstantiated,
//                            DeltaTests,
//                            testing::ValuesIn(GetTestParameters()));
//}