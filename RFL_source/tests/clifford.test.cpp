//
// Created by Paul Druce on 11/08/2022.
//
#include "Cliff.hpp"
#include <gtest/gtest.h>

TEST(CliffordTest, NoErrorWhenConstructing) {
  typedef struct {
    int p;
    int q;
  } CliffData;

  std::vector<CliffData> cliff_data = {
      {1, 0},
      {0, 1},
      {1, 1},
      {1, 2},
      {2, 1},
      {1, 3},
      {3, 1},
      {3, 2},
      {2, 3},
      {3, 3}};

  for (auto& d : cliff_data) {
    EXPECT_NO_THROW(
        Cliff cliff(d.p, d.q););
  }
}

TEST(CliffordTest, GetP) {
  for (int p = 1; p < 8; p++) {
    Cliff cliff(p, 1);

    EXPECT_EQ(p, cliff.get_p());
  }
}

TEST(CliffordTest, GetQ) {
  for (int q = 1; q < 8; q++) {
    Cliff cliff(1, q);

    EXPECT_EQ(q, cliff.get_q());
  }
}

TEST(CliffordTest, GetDimGamma) {

  typedef struct {
    int p;
    int q;
    int dim_gamma;
  } DimGammaData;

  std::vector<DimGammaData> data = {
      {1, 0, 1},
      {0, 1, 1},
      {1, 1, 2},
      {1, 2, 2},
      {2, 1, 2},
      {1, 3, 4},
      {3, 1, 4},
      {3, 3, 8}};
  for (auto& d : data) {
    Cliff cliff(d.p, d.q);
    EXPECT_EQ(d.dim_gamma, cliff.get_dim_gamma()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(CliffordTest, GetGammas) {
  typedef struct {
    int p;
    int q;
    int num_gammas;
  } GetGammasData;

  std::vector<GetGammasData> data = {
      {1, 0, 1},
      {0, 1, 1},
      {1, 1, 2},
      {1, 2, 3},
      {2, 1, 3},
      {1, 3, 4},
      {3, 3, 6}};

  for (auto& d : data) {
    Cliff cliff(d.p, d.q);
    EXPECT_EQ(d.num_gammas, cliff.get_gamma().size()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}
