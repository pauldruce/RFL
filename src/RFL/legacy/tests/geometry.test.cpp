//
// Created by Paul Druce on 11/08/2022.
//
#include "Geom24.hpp"
#include <gtest/gtest.h>

TEST(GeometryTest, NoErrorWhenConstructing) {
  typedef struct GeomData {
    int p;
    int q;
    int dim;
    double g2;
  } GeomData;

  std::vector<GeomData> geom_data = {
      {1, 1, 5, -2.7},
      {1, 3, 6, -3.54},
      {2, 1, 7, -1.57}};

  for (auto& data : geom_data) {
    EXPECT_NO_THROW(
        Geom24 geom(data.p, data.q, data.dim, data.g2));
  }
}

TEST(GeometryTest, GetP) {
  for (int p = 1; p < 5; p++) {
    Geom24 geom(p, 1, 5, -2.7);

    EXPECT_EQ(p, geom.get_p());
  }
}

TEST(GeometryTest, GetQ) {
  for (int q = 1; q < 5; q++) {
    Geom24 geom(1, q, 5, -2.5);

    EXPECT_EQ(q, geom.get_q());
  }
}

TEST(GeometryTest, GetDim) {
  for (int d = 2; d < 10; d++) {
    Geom24 geom(2, 3, d, -2.7);

    EXPECT_EQ(d, geom.get_dim());
  }
}

TEST(GeometryTest, GetNumHMatrices) {
  typedef struct {
    int p;
    int q;
    int num_h_mat;
  } NumHMat_data;

  std::vector<NumHMat_data> data = {
      {1, 1, 1},
      {1, 2, 1},
      {2, 1, 3},
      {3, 3, 16}};

  for (auto& d : data) {
    Geom24 geom(d.p, d.q, 5, -2.7);
    EXPECT_EQ(d.num_h_mat, geom.get_nH()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(GeometryTest, GetNumLMatrices) {
  typedef struct {
    int p;
    int q;
    int num_L_mat;
  } NumLMat_data;

  std::vector<NumLMat_data> data = {
      {1, 1, 1},
      {1, 2, 3},
      {2, 1, 1},
      {3, 3, 16}};

  for (auto& d : data) {
    Geom24 geom(d.p, d.q, 5, -2.7);
    EXPECT_EQ(d.num_L_mat, geom.get_nL()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(GeometryTest, GetNumHLMatrices) {
  typedef struct {
    int p;
    int q;
    int num_hl_mat;
  } NumHLMat_data;

  std::vector<NumHLMat_data> data = {
      {1, 1, 2},
      {1, 2, 4},
      {2, 1, 4},
      {3, 3, 32}};

  for (auto& d : data) {
    Geom24 geom(d.p, d.q, 5, -2.7);
    EXPECT_EQ(d.num_hl_mat, geom.get_nHL()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}
