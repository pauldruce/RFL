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

  const std::vector<NumHMatData> data = {
      {1, 1, 1},
      {1, 2, 1},
      {2, 1, 3},
      {3, 3, 16}};

  for (const auto& [p, q, num_h_mat] : data) {
    DiracOperator dirac(p, q, 5);
    EXPECT_EQ(num_h_mat, dirac.getNumHermitianMatrices()) << "(p,q) = (" << p << "," << q << ")";
  }
}

TEST(DiracOperatorTests, NumAntiHermitianMatricesIsCorrectlySet) {
  typedef struct {
    int p;
    int q;
    int num_h_mat;
  } NumLMatData;

  const std::vector<NumLMatData> data = {
      {1, 1, 1},
      {1, 2, 3},
      {2, 1, 1},
      {3, 3, 16}};

  for (const auto& [p, q, num_h_mat] : data) {
    DiracOperator dirac(p, q, 5);
    EXPECT_EQ(num_h_mat, dirac.getNumAntiHermitianMatrices()) << "(p,q) = (" << p << "," << q << ")";
  }
}

TEST(DiracOperatorTests, TotalNumMatricesIsSetCorrectly) {
  typedef struct {
    int p;
    int q;
    int num_hl_mat;
  } NumTotalMatricesData;

  const std::vector<NumTotalMatricesData> data = {
      {1, 1, 2},
      {1, 2, 4},
      {2, 1, 4},
      {3, 3, 32}};

  for (const auto& [p, q, num_hl_mat] : data) {
    DiracOperator dirac(p, q, 5);
    EXPECT_EQ(num_hl_mat, dirac.getNumMatrices()) << "(p,q) = (" << p << "," << q << ")";
  }
}

TEST(DiracOperatorTests, GetType) {
  const std::vector<std::pair<int, int>> data = {
      {1, 1},
      {1, 2},
      {2, 1},
      {3, 3}};

  for (const auto& [p, q] : data) {
    DiracOperator dirac(p, q, 5);
    EXPECT_EQ(p, dirac.getType().first);
    EXPECT_EQ(q, dirac.getType().second);
  }
}

TEST(DiracOperatorTests, GetEigenvalues) {
  const std::vector<std::pair<int, int>> data = {
      {1, 1},
      {1, 2},
      {2, 1},
      {3, 3}};

  for (const auto& [p, q] : data) {
    DiracOperator dirac(p, q, 5);

    auto eigenvalues = dirac.getEigenvalues();

    auto dirac_matrix = dirac.getDiracMatrix();
    auto expected_eigenvalues = arma::eig_sym(dirac_matrix);

    ASSERT_EQ(eigenvalues.n_elem, expected_eigenvalues.n_elem);

    eigenvalues = arma::sort(eigenvalues);
    expected_eigenvalues = arma::sort(expected_eigenvalues);

    for (unsigned i = 0; i < (unsigned)eigenvalues.n_elem; i++) {
      EXPECT_FLOAT_EQ(eigenvalues[i], expected_eigenvalues[i]);
    }
  }
}

TEST(DiracOperatorTests, CopyConstructorWorks) {
  const DiracOperator dirac_1(1, 3, 5);
  const DiracOperator dirac_2 = dirac_1;// NOLINT(*-unnecessary-copy-initialization)

  ASSERT_TRUE(arma::approx_equal(
      dirac_1.getDiracMatrix(),
      dirac_2.getDiracMatrix(),
      "absdiff",
      1e-10));
}

TEST(DiracOperatorTests, GetHermitianMatricesReturnsHermitianMatrices) {
  typedef struct {
    int p;
    int q;
    int dim;
    int n_herm_matrices;
  } DiracParams;

  const std::vector<DiracParams> params{
      {2, 1, 5, 3},
      {1, 2, 6, 1},
      {3, 3, 7, 16}};

  for (const auto& [p, q, dim, n_herm_matrices] : params) {
    DiracOperator dirac(p, q, dim);

    auto hermitian_matrices = dirac.getHermitianMatrices();
    EXPECT_EQ(n_herm_matrices, hermitian_matrices.size())
        << "For dirac params"
        << "(p,q)=(" << p << "," << q << ")"
        << " dim = " << dim;

    for (auto& matrix : hermitian_matrices) {
      EXPECT_TRUE(matrix.is_hermitian())
          << "For dirac params"
          << "(p,q)=(" << p << "," << q << ")"
          << " dim = " << dim;
    }
  }
}

TEST(DiracOperatorTests, GetAntiHermitianMatricesReturnsAntiHermitianMatrices) {
  typedef struct {
    int p;
    int q;
    int dim;
    int n_anti_herm_matrices;
  } DiracParams;

  const std::vector<DiracParams> params{
      {2, 1, 5, 1},
      {1, 2, 6, 3},
      {3, 3, 7, 16}};

  for (const auto& [p, q, dim, n_anti_herm_matrices] : params) {
    DiracOperator dirac(p, q, dim);

    auto hermitian_matrices = dirac.getAntiHermitianMatrices();
    EXPECT_EQ(n_anti_herm_matrices, hermitian_matrices.size())
        << "For dirac params"
        << "(p,q)=(" << p << "," << q << ")"
        << " dim = " << dim;

    for (auto& matrix : hermitian_matrices) {
      EXPECT_TRUE(!matrix.is_hermitian())
          << "For dirac params"
          << "(p,q)=(" << p << "," << q << ")"
          << " dim = " << dim;
    }
  }
}