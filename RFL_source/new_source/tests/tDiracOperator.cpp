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

TEST(DiracOperatorTests,GetMoms){
  DiracOperator D(1,1,5);
  arma::cx_mat identityMatrix(5,5);
  auto* moms = D.get_moms();
  identityMatrix.eye();

  for(int i = 0; i<D.nHL; i++) {
	ASSERT_TRUE(
		arma::approx_equal(moms[i], identityMatrix, "absdiff", 1e-10)
	);
  }

}

//TEST(DiracOperatorTests, BuildDiracHasCorrectDims) {
//  std::vector<DiracOpData> data = {
//	  {1, 0, 5},
//	  {0, 1, 5},
//	  {1, 3, 6},
//	  {3, 1, 7},
//	  {2, 2, 10}
//  };
//
//  for (auto &d : data) {
//
//  }
//}