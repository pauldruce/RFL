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
  std::vector<CliffordData> cliff_data = {
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
        Clifford cliff(d.p, d.q););
  }
}

TEST(CliffordTests, GetP) {
  for (int p = 1; p < 8; p++) {
    Clifford C(p, 1);

    EXPECT_EQ(p, C.get_p());
  }
}

TEST(CliffordTests, GetQ) {
  for (int q = 1; q < 8; q++) {
    Clifford C(1, q);
    EXPECT_EQ(q, C.get_q());
  }
}

TEST(CliffordTests, GetDimGamma) {
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
    Clifford cliff(d.p, d.q);
    EXPECT_EQ(d.dim_gamma, cliff.get_dim_gamma()) << "(p,q) = (" << d.p << "," << d.q << ")";
  }
}

TEST(CliffordTests, GetGammas) {
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
    Clifford cliff(d.p, d.q);
    EXPECT_EQ(d.num_gammas, cliff.get_gammas().size()) << "(p,q) = (" << d.p << "," << d.q << ")";
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

  for (const auto& d : data) {
    Clifford C(d.p, d.q);

    auto gammas = C.get_gammas();
    std::vector<arma::cx_mat> herm_gammas(gammas.begin(), gammas.begin() + d.p);
    std::vector<arma::cx_mat> anti_herm_gammas(gammas.begin() + d.p, gammas.end());

    EXPECT_EQ(d.p, herm_gammas.size());
    for (auto& hg : herm_gammas) {
      EXPECT_TRUE(hg.is_hermitian()) << "(p,q) = (" << d.p << "," << d.q << ")\n ";
    }

    EXPECT_EQ(d.q, anti_herm_gammas.size());
    for (auto& ahg : anti_herm_gammas) {
      EXPECT_TRUE(!ahg.is_hermitian()) << "(p,q) = (" << d.p << "," << d.q << ")";
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

  for (const auto& d : data) {
    Clifford C(d.p, d.q);
    int n = d.p + d.q;
    int exponent = (n % 2 == 0) ? (int)n / 2 : (int)((n - 1) / 2);
    int expect_dim = 1 << exponent;

    // Easy check
    EXPECT_EQ(expect_dim, C.get_dim_gamma());

    // Explicit check for all gammas
    auto gammas = C.get_gammas();
    for (auto& g : gammas) {
      ASSERT_EQ(g.n_rows, g.n_cols)
                  << "Gamma matrices not squarefailed for (p,q) = (" << d.p << "," << d.q << ")";
      EXPECT_EQ(expect_dim, g.n_rows) << "Gamma dim not correct for (p,q) = (" << d.p << "," << d.q << ")";
    }
  }
}

// TODO: find minimum set of configurations of (p,q) that covers all code.
TEST(CliffordTests, ChiralityIsCorrect) {
  // TODO: Fix source code, because this test fails if max_p, max_q is >= 6
  constexpr int max_p = 5;
  constexpr int max_q = 5;
  vector<CliffordData> data = {};

  data.push_back({1, 0});
  data.push_back({0, 1});
  for (int q = 1; q < max_p; q++) {
    for (int p = 1; p < max_q; p++) {
      data.push_back({p, q});
    }
  }

  for (auto& d : data) {
    Clifford C(d.p, d.q);
    int s = (d.q - d.p + 8 * d.p) % 8;// Need to add 8*p because % does behave for negative values.
    auto gammas = C.get_gammas();
    auto chirality = C.get_chiral();

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
  Clifford C2(2, 1);

  C1 *= C2;

  EXPECT_EQ(C1.get_p(), 3);
  EXPECT_EQ(C1.get_q(), 3);
}

TEST(CliffordTests, CliffordModuleMultiplication) {
  Clifford C1(1, 2);
  Clifford C2(2, 1);

  auto C3 = C1 * C2;

  // New Clifford module has correct dimensions
  EXPECT_EQ(C3.get_p(), 3);
  EXPECT_EQ(C3.get_q(), 3);

  // Originals are not modified
  EXPECT_EQ(C1.get_p(), 1);
  EXPECT_EQ(C1.get_q(), 2);

  EXPECT_EQ(C2.get_p(), 2);
  EXPECT_EQ(C2.get_q(), 1);
}
// TODO: Write a test to check the first few types have correct gammas
