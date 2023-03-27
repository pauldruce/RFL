//
// Created by Paul Druce on 28/12/2022.
//

#include "DiracOperator.hpp"
#include <gtest/gtest.h>

TEST(DiracOperatorTests, NoErrorsWhenConstruction) {
  for (int p = 0; p < 5; p++) {
    for (int q = 0; q < 5; q++) {
      for (int dim = 5; dim < 8; dim++) {
        if (p + q > 5 || (p == 0 && q == 0)) {
          continue;
        }
        EXPECT_NO_THROW(DiracOperator dirac(p, q, dim)) << "(p,q,dim) = (" << p << "," << q << "," << dim << ")";
      }
    }
  }
}

TEST(DiracOperatorTests, NumHermitianMatricesIsCorrectlySet) {
  typedef struct {
    int p;
    int q;
    int num_h_mat;
  } NumHMatData;

  std::vector<NumHMatData> data = {
      {1, 1, 1},
      {1, 2, 1},
      {2, 1, 3},
      {3, 3, 16}};

  for (auto& d : data) {
    DiracOperator dirac(d.p, d.q, 5);
    EXPECT_EQ(d.num_h_mat, dirac.getNumHermitianMatrices()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(DiracOperatorTests, NumAntiHermitianMatricesIsCorrectlySet) {
  typedef struct {
    int p;
    int q;
    int num_h_mat;
  } NumLMatData;

  std::vector<NumLMatData> data = {
      {1, 1, 1},
      {1, 2, 3},
      {2, 1, 1},
      {3, 3, 16}};

  for (auto& d : data) {
    DiracOperator dirac(d.p, d.q, 5);
    EXPECT_EQ(d.num_h_mat, dirac.getNumAntiHermitianMatrices()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(DiracOperatorTests, TotalNumMatricesIsSetCorrectly) {
  typedef struct {
    int p;
    int q;
    int num_hl_mat;
  } NumTotalMatricesData;

  std::vector<NumTotalMatricesData> data = {
      {1, 1, 2},
      {1, 2, 4},
      {2, 1, 4},
      {3, 3, 32}};

  for (auto& d : data) {
    DiracOperator dirac(d.p, d.q, 5);
    EXPECT_EQ(d.num_hl_mat, dirac.getNumMatrices()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}