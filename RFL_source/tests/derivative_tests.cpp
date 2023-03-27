//
// Created by Paul Druce on 01/07/2022.
//
#include "Geom24.hpp"
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>
#include <iostream>

namespace DerivativeTests {
using namespace std;
using namespace arma;

TEST(DerivativeTests, DerivativesAreWhatExactly) {
  typedef struct DerivativeParams {
    int p;
    int q;
    int dim;
    double g2;
  } DerivativeParams;

  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));

  // Create (p,q) combo's upto n = p+q=5.
  std::pair<int, int> pq_pairs[] = {
      {1, 0},
      {0, 1},
      {2, 0},
      {1, 1},
      {0, 2}};

  for (auto pair : pq_pairs) {
    for (int dim = 2; dim < 8; dim++) {
      for (double g2 = -5.0; g2 < 0.0; g2 += 1.0) {
        DerivativeParams params = {pair.first, pair.second, dim, g2};

        Geom24 G(params.p, params.q, params.dim, params.g2);
        ifstream in("derivative_date.txt");
        G.read_mat(in);
        for (int i = 0; i < G.get_nHL(); ++i) {
          EXPECT_TRUE(G.get_mat(i).is_diagmat()) << "Matrix is not diagonal i=" << i << ", " << G.get_mat(i);
        }

        EXPECT_TRUE(G.der_dirac4(0, true).is_diagmat()) << "G.derDirac4(0,true) matrix is not diagonal";
      }
    }
  }

  gsl_rng_free(engine);
}

}// namespace DerivativeTests
 //
 //   struct DerivativeParams{
 //      int p;
 //      int q;
 //      int dim;
 //      double g2;
 //   };
 //
 //   class DerivativeTests : public testing::TestWithParam<DerivativeParams>
 //   {
 //      void SetUp() override
 //      {
 //         gsl_rng *engine = gsl_rng_alloc(gsl_rng_ranlxd1);
 //         gsl_rng_set(engine, time(NULL));
 //      }
 //   };
 //
 //   TEST_P(DerivativeTests, DerivativesAreWhatsExactly)
 //   {
 //      DerivativeParams params = GetParam();
 //      Geom24 G(params.p, params.q, params.dim, params.g2);
 //      ifstream in("derivative_date.txt");
 //      G.read_mat(in);
 //      for (int i = 0; i < G.get_nHL(); ++i)
 //         cout << G.get_mat(i) << endl;
 //
 //      cout << G.der_dirac4(0, true) << endl;
 //
 //      // As I'm not sure what these Tests should be checking.
 //      // I'm making them fail until I've spoken with Mauro.
 //      ASSERT_TRUE(false);
 //   }
 //
 //   std::vector<DerivativeParams> GetTestParams()
 //   {
 //      std::vector<DerivativeParams> test_params;
 //
 //      // Create (p,q) combo's upto n = p+q=5.
 //      std::pair<int, int> pq_pairs[] = {
 //            {1, 0},{0, 1},
 //            {2, 0},{1, 1},{0, 2}
 //      };
 //
 //      for (auto pair: pq_pairs)
 //      {
 //         for (int dim = 2; dim < 10; dim++)
 //         {
 //            for (double g2 = -5.0; g2 < 0.0; g2 += 1.0)
 //            {
 //               DerivativeParams params = {pair.first, pair.second, dim, g2};
 //               test_params.emplace_back(params);
 //            }
 //         }
 //      }
 //      return test_params;
 //   }
 //
 //   INSTANTIATE_TEST_SUITE_P(DerivativeTestsInstantiated,
 //                            DerivativeTests,
 //                            testing::ValuesIn(GetTestParams()));
 //}