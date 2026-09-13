//
// Created by Paul Druce on 28/06/2022.
//

#include "Geom24.hpp"
#include <gtest/gtest.h>
namespace CompilationTests {
TEST(CompilationTests, Geom24CompilesWithoutError) {
  struct CompilationParameters {
    int p;
    int q;
    int dim;
    double g2;
  };

  std::vector<CompilationParameters> test_params;

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
    for (int dim = 2; dim < 8; dim++) {
      for (double g2 = -5.0; g2 < 0.0; g2 += 1.0) {
        CompilationParameters params = {pair.first, pair.second, dim, g2};
        test_params.emplace_back(params);
      }
    }
  }

  for (auto& param : test_params) {
    ASSERT_NO_THROW(Geom24(param.p, param.q, param.dim, param.g2)) << "";
  }
}

// PARAMETRISED TEST USING GTEST
//
//   struct CompilationParameters
//   {
//      int p;
//      int q;
//      int dim;
//      double g2;
//   };
//   class CompilationTests : public ::testing::TestWithParam<CompilationParameters>
//   {
//   };
//
//   TEST_P(CompilationTests, Geom24CompilesWithoutError)
//   {
//      CompilationParameters params = GetParam();
//      ASSERT_NO_THROW(Geom24(params.p, params.q, params.dim, params.g2)) << "";
//   }
//
//   static std::vector<CompilationParameters> GetTestParameters()
//   {
//      std::vector<CompilationParameters> test_params;
//
//      // Create (p,q) combo's upto n = p+q=5.
//      std::pair<int, int> pq_pairs[] = {
//            {1, 0},
//            {0, 1},
//            {2, 0},
//            {1, 1},
//            {0, 2},
//            {3, 0},
//            {2, 1},
//            {1, 2},
//            {0, 3},
//            {4, 0},
//            {3, 1},
//            {2, 2},
//            {1, 3},
//            {0, 4},
//            {5, 0},
//            {4, 1},
//            {3, 2},
//            {2, 3},
//            {1, 4},
//            {0, 5}
//      };
//
//      for (auto pair: pq_pairs)
//      {
//         for (int dim = 2; dim < 10; dim++)
//         {
//            for (double g2 = -5.0; g2 < 0.0; g2 += 1.0)
//            {
//               CompilationParameters params = {pair.first, pair.second, dim, g2};
//               test_params.emplace_back(params);
//            }
//         }
//      }
//      return test_params;
//   }
//
//   INSTANTIATE_TEST_SUITE_P(FastCompilationTests,
//                            CompilationTests,
//                            testing::ValuesIn(GetTestParameters()));
}// namespace CompilationTests
