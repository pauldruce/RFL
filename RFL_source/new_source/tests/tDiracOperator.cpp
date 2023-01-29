//
// Created by Paul Druce on 28/12/2022.
//


#include <gtest/gtest.h>
#include "DiracOperator.hpp"

typedef struct {
  int p;
  int q;
  int dim;
} DiracOpData;

TEST(DiracOperatorTests, NoErrorsWhenConstruction) {
  for (int p = 0; p < 5; p++) {
	for (int q = 0; q < 5; q++) {
	  for (int dim = 5; dim < 8; dim++) {
		if (p + q > 5 || (p == 0 && q == 0)) {
		  continue;
		}
		EXPECT_NO_THROW(DiracOperator D(p, q, dim)) << "(p,q,dim) = (" << p << "," << q << "," << dim << ")";
	  }
	}
  }
}

TEST(DiracOperatorTests,NumH_IsCorrectlySet){
  typedef struct
  {
	int p;
	int q;
	int num_h_mat;
  } NumHMat_data;

  std::vector<NumHMat_data> data = {
	  {1, 1, 1},
	  {1, 2, 1},
	  {2, 1, 3},
	  {3, 3, 16}};

  for (auto &d : data)
  {
	DiracOperator D(d.p, d.q, 5);
	EXPECT_EQ(d.num_h_mat, D.nH) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(DiracOperatorTests,NumL_IsCorrectlySet){
  typedef struct
  {
	int p;
	int q;
	int num_h_mat;
  } NumLMat_data;

  std::vector<NumLMat_data> data = {
	  {1, 1, 1},
	  {1, 2, 3},
	  {2, 1, 1},
	  {3, 3, 16}};

  for (auto &d : data)
  {
	DiracOperator D(d.p, d.q, 5);
	EXPECT_EQ(d.num_h_mat, D.nL) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(DiracOperatorTests,NumHL_IsCorrectlySet){
  typedef struct
  {
	int p;
	int q;
	int num_hl_mat;
  } NumHLMat_data;

  std::vector<NumHLMat_data> data = {
	  {1, 1, 2},
	  {1, 2, 4},
	  {2, 1, 4},
	  {3, 3, 32}};

  for (auto &d : data)
  {
	DiracOperator D(d.p, d.q, 5);
	EXPECT_EQ(d.num_hl_mat, D.nHL) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}