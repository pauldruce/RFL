//
// Created by Paul Druce on 27/12/2022.
//
#include <gtest/gtest.h>
#include "Action.hpp"

static auto engine = gsl_rng_alloc(gsl_rng_ranlxd1);

static void CompareActions(int p, int q, int dim, double g2) {
  constexpr int numOfTestRepeats = 100;
  Action A(g2);
  DiracOperator D(p, q, dim);

  for (int i = 0; i < numOfTestRepeats; ++i) {
	D.randomise(engine);
	double d2 = D.dim * D.dim;

	auto S1 = A.calculate_S(D) / d2;

	auto S2 = A.calculate_S_from_dirac(D) / d2;

	EXPECT_TRUE(fabs(S1 - S2) < 1e-8) << "Methods differ more then 1e-8";
  }
}

TEST(ActionTests, ActionMethodsDontDiffer) {
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

  for (auto &d : testing_params) {
	CompareActions(d.p, d.q, d.dim, d.g2);
  }
}

TEST(ActionTests, Set_g2) {
  const double initial_g2 = 1.1;
  Action A(initial_g2);
  EXPECT_EQ(A.get_g2(), initial_g2);

  const double new_g2 = -3.4;
  A.set_g2(new_g2);
  EXPECT_EQ(A.get_g2(), new_g2);
}

TEST(ActionTests, Set_g4) {
  const double initial_g4 = 1.1;
  Action A(1.0, initial_g4);
  EXPECT_EQ(A.get_g4(), initial_g4);

  const double new_g4 = -3.4;
  A.set_g4(new_g4);
  EXPECT_EQ(A.get_g4(), new_g4);
}

TEST(ActionTests, SetParameters) {
  const double initial_g2 = 1.1, initial_g4 = 2.2;
  Action A(initial_g2, initial_g4);

  EXPECT_EQ(A.get_g2(), initial_g2);
  EXPECT_EQ(A.get_g4(), initial_g4);

  const double new_g2 = -5.0, new_g4 = 7.2;
  A.set_params(new_g2, new_g4);

  EXPECT_EQ(A.get_g2(), new_g2);
  EXPECT_EQ(A.get_g4(), new_g4);
}

TEST(ActionTests, CreateWithNoParams){
  EXPECT_NO_THROW(
	  Action A;
	  EXPECT_EQ(A.get_g2(), 0);
	  EXPECT_EQ(A.get_g4(), 0);
  );
}

TEST(ActionTests, PrintAction){
 DiracOperator D(1,1,5);
 Action A(2.0,4.0);

 std::stringstream capturedStream;
 A.print_S(D,capturedStream);

 EXPECT_EQ(capturedStream.str(), "200 800\n");
}