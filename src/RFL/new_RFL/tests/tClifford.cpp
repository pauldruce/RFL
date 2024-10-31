//
// Created by Paul Druce on 27/12/2022.
//

#include "Clifford.hpp"
#include <gtest/gtest.h>

using namespace arma;
using namespace std;

typedef struct {
  int p;
  int q;
} CliffordData;

TEST(CliffordTests, NoErrorWhenConstructing) {
  const std::vector<CliffordData> cliff_data = {
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

  for (const auto& [p, q] : cliff_data) {
    EXPECT_NO_THROW(
        Clifford cliff(p, q););
  }
}

TEST(CliffordTests, GetP) {
  for (int p = 1; p < 8; p++) {
    Clifford clifford(p, 1);

    EXPECT_EQ(p, clifford.getP());
  }
}

TEST(CliffordTests, GetQ) {
  for (int q = 1; q < 8; q++) {
    Clifford clifford(1, q);
    EXPECT_EQ(q, clifford.getQ());
  }
}

TEST(CliffordTests, GetDimGamma) {
  typedef struct {
    int p;
    int q;
    int dim_gamma;
  } DimGammaData;

  const std::vector<DimGammaData> data = {
      {1, 0, 1},
      {0, 1, 1},
      {1, 1, 2},
      {1, 2, 2},
      {2, 1, 2},
      {1, 3, 4},
      {3, 1, 4},
      {3, 3, 8}};
  for (const auto& [p, q, dim_gamma] : data) {
    Clifford cliff(p, q);
    EXPECT_EQ(dim_gamma, cliff.getGammaDimension()) << "(p,q) = (" << p << "," << q << ")";
  }
}

TEST(CliffordTests, GetGammas) {
  typedef struct {
    int p;
    int q;
    int num_gammas;
  } GetGammasData;

  const std::vector<GetGammasData> data = {
      {1, 0, 1},
      {0, 1, 1},
      {1, 1, 2},
      {1, 2, 3},
      {2, 1, 3},
      {1, 3, 4},
      {3, 3, 6}};

  for (const auto& [p, q, num_gammas] : data) {
    Clifford cliff(p, q);
    EXPECT_EQ(num_gammas, cliff.getGammaMatrices().size()) << "(p,q) = (" << p << "," << q << ")";
  }
}

TEST(CliffordTests, GammasHaveCorrectHermiticity) {
  vector<CliffordData> data;
  data.push_back({1, 0});
  data.push_back({0, 1});
  for (int p = 1; p < 8; p++) {
    for (int q = 1; q < 8; q++) {
      data.push_back({p, q});
    }
  }

  for (const auto& [p, q] : data) {
    Clifford C(p, q);

    auto gammas = C.getGammaMatrices();
    std::vector herm_gammas(gammas.begin(), gammas.begin() + p);
    std::vector anti_herm_gammas(gammas.begin() + p, gammas.end());

    EXPECT_EQ(p, herm_gammas.size());
    for (auto& hg : herm_gammas) {
      EXPECT_TRUE(hg.is_hermitian()) << "(p,q) = (" << p << "," << q << ")\n ";
    }

    EXPECT_EQ(q, anti_herm_gammas.size());
    for (auto& ahg : anti_herm_gammas) {
      EXPECT_TRUE(!ahg.is_hermitian()) << "(p,q) = (" << p << "," << q << ")";
    }
  }
}

TEST(CliffordTests, GammasHaveCorrectDims) {

  vector<CliffordData> data;
  data.push_back({1, 0});
  data.push_back({0, 1});
  for (int p = 1; p < 8; p++) {
    for (int q = 1; q < 8; q++) {
      data.push_back({p, q});
    }
  }

  for (const auto& [p, q] : data) {
    Clifford C(p, q);
    const int n = p + q;
    const int exponent = (n % 2 == 0) ? (int)n / 2 : (int)((n - 1) / 2);
    const int expect_dim = 1 << exponent;

    // Easy check
    EXPECT_EQ(expect_dim, C.getGammaDimension());

    // Explicit check for all gammas
    auto gammas = C.getGammaMatrices();
    for (auto& g : gammas) {
      ASSERT_EQ(g.n_rows, g.n_cols)
          << "Gamma matrices not squarefailed for (p,q) = (" << p << "," << q << ")";
      EXPECT_EQ(expect_dim, g.n_rows) << "Gamma dim not correct for (p,q) = (" << p << "," << q << ")";
    }
  }
}

// TODO: find minimum set of configurations of (p,q) that covers all code.
TEST(CliffordTests, ChiralityIsCorrect) {
  // TODO: Fix source code, because this test fails if max_p, max_q is >= 6
  constexpr int max_p = 5;
  vector<CliffordData> data = {};

  data.push_back({1, 0});
  data.push_back({0, 1});
  for (int q = 1; q < max_p; q++) {
    constexpr int max_q = 5;
    for (int p = 1; p < max_q; p++) {
      data.push_back({p, q});
    }
  }

  for (const auto& [p, q] : data) {
    Clifford C(p, q);
    const int s = (q - p + 8 * p) % 8;// Need to add 8*p because % does behave for negative values.
    auto gammas = C.getGammaMatrices();
    auto chirality = C.getChiral();

    int exponent = (int)s * (s + 1) / 2;
    arma::cx_mat gammas_producted = arma::eye<arma::cx_mat>(gammas[0].n_rows, gammas[0].n_cols);
    for (auto& g : gammas) {
      gammas_producted = gammas_producted * g;
    }
    arma::cx_mat
        expect_chirality =
            gammas_producted * (std::complex<double>)std::pow(std::complex<double>(0, 1), exponent);
    expect_chirality.clean(1e-10);
    EXPECT_TRUE(arma::approx_equal(expect_chirality, chirality, "absdiff", 1e-4))
        << "Chirality Operator not correct for " << C
        << "\nExpected Chirality:\n"
        << expect_chirality
        << "\nRetrieved Chirality:\n"
        << chirality;
  }
}

TEST(CliffordTests, CliffordModuleInlineMultiplication) {
  Clifford C1(1, 2);
  const Clifford C2(2, 1);

  C1 *= C2;

  EXPECT_EQ(C1.getP(), 3);
  EXPECT_EQ(C1.getQ(), 3);
}

TEST(CliffordTests, CliffordModuleMultiplication) {
  const Clifford C1(1, 2);
  const Clifford C2(2, 1);

  const auto C3 = C1 * C2;

  // New Clifford module has correct dimensions
  EXPECT_EQ(C3.getP(), 3);
  EXPECT_EQ(C3.getQ(), 3);

  // Originals are not modified
  EXPECT_EQ(C1.getP(), 1);
  EXPECT_EQ(C1.getQ(), 2);

  EXPECT_EQ(C2.getP(), 2);
  EXPECT_EQ(C2.getQ(), 1);
}
// TODO: Write a test to check the first few types have correct gammas
