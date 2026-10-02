//
// Unit tests verifying analytic trace variations against direct action differences.
//

#include "../BarrettGlaser/BarrettGlaserAction.hpp"
#include "../DiracOperator.hpp"
#include "../StdRng.hpp"
#include <armadillo>
#include <complex>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

namespace {

// Absolute numerical tolerance for double-precision trace differences.
// Accounts for O(N^3) matrix trace cancellations and condition numbers;
// kTolerance = 1e-7 safely bounds floating-point cancellation across random matrix perturbations.
constexpr double kTolerance = 1e-7;

const std::vector<std::pair<int, int>> kGeometries = {
    {1, 0},
    {0, 1},
    {1, 1},
    {2, 0},
    {0, 2},
    {2, 2},
    {1, 3},
    {4, 0}};

// Applies an elementary variation delta_M to matrix x at entry (row, col).
//
// 1. Invariant:
//    Hermiticity preservation: (M + delta_M)† == M + delta_M.
//
// 2. Mathematical Move:
//    The Barrett-Glaser analytic trace variations evaluate the change under:
//    delta_M = z * E_{row,col} + conj(z) * E_{col,row}
//    - Off-diagonal (row != col): updates both (row, col) and (col, row).
//    - Diagonal (row == col): collapses to z + conj(z) = 2 * Re(z). The imaginary part
//      cancels identically for any z in C, ensuring the diagonal variation is strictly real.
//
// 3. Dual-Path Verification:
//    Used to perturb a copied Dirac operator to verify analytic delta formulas.
void applyPerturbation(
    DiracOperator& dirac,
    const int x,
    const int row,
    const int col,
    const arma::cx_double& z) {
  auto& mat = dirac.getMatrix(x);
  if (row != col) {
    mat(row, col) += z;
    mat(col, row) += std::conj(z);
  } else {
    // For diagonal entries, delta_M = z + conj(z) = 2 * Re(z).
    // Using 2.0 * z.real() ensures the diagonal remains strictly real even if
    // an imaginary component were passed.
    mat(row, row) += 2.0 * z.real();
  }
}

// Generates an elementary perturbation value z for matrix entry (row, col).
//
// 1. Invariant:
//    Diagonal entries of Hermitian matrices must remain strictly real (M_{ii} == conj(M_{ii})).
//    Setting Im(z) = 0 prevents the perturbation from violating Hermiticity.
//
// 2. Off-Diagonal Entries:
//    Can take arbitrary complex values because applyPerturbation() updates the
//    transposed symmetric entry with conj(z), maintaining Hermiticity.
arma::cx_double randomPerturbation(
    const IRng& rng,
    const int row,
    const int col) {
  const double re = -0.5 + rng.getUniform();
  if (row == col) {
    // Diagonal entries of Hermitian matrices must remain strictly real.
    return {re, 0.0};
  }
  const double im = -0.5 + rng.getUniform();
  return {re, im};
}

// Computes the direct numerical difference:
//   ΔTr(D²) = Tr((D')²) - Tr(D²)
//
// Dual-path verification ground truth for BarrettGlaserAction::delta2.
double calculateDirectDelta2(
    DiracOperator dirac,
    const int x,
    const int row,
    const int col,
    const arma::cx_double& z) {
  const double initial_trace = dirac.traceOfDiracSquared();
  applyPerturbation(dirac, x, row, col, z);
  const double final_trace = dirac.traceOfDiracSquared();
  return final_trace - initial_trace;
}

// Computes the direct numerical difference:
//   ΔTr(D⁴) = Tr((D')⁴) - Tr(D⁴)
//
// Dual-path verification ground truth for BarrettGlaserAction::delta4.
double calculateDirectDelta4(
    DiracOperator dirac,
    const int x,
    const int row,
    const int col,
    const arma::cx_double& z) {
  const double initial_trace = dirac.traceOfDirac4();
  applyPerturbation(dirac, x, row, col, z);
  const double final_trace = dirac.traceOfDirac4();
  return final_trace - initial_trace;
}

// Computes the direct numerical action difference:
//   ΔS = S(D') - S(D) = g_2 * ΔTr(D²) + g_4 * ΔTr(D⁴)
//
// Dual-path verification ground truth for BarrettGlaserAction::delta24.
double calculateDirectDelta24(
    const BarrettGlaserAction& action,
    DiracOperator dirac,
    const int x,
    const int row,
    const int col,
    const arma::cx_double& z) {
  const double initial_action = action.calculateS(dirac);
  applyPerturbation(dirac, x, row, col, z);
  const double final_action = action.calculateS(dirac);
  return final_action - initial_action;
}

}// namespace

// =============================================================================
// Unit Tests
// =============================================================================

// 1. Invariant: Action variation under Hermitian rank-1 update Tr(D^2(M + delta_M)) - Tr(D^2(M)).
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Eq. (3.1).
// 3. Dual-Path Verification: BarrettGlaserAction::delta2 vs calculateDirectDelta2.
// 4. Tolerance Rationale: kTolerance (1e-7) bounds double-precision trace cancellation.
TEST(DeltaTests, Delta2MatchesTraceDifference) {
  const StdRng rng(42);

  for (const auto& [p, q] : kGeometries) {
    for (int dim = 2; dim <= 4; ++dim) {
      DiracOperator dirac(p, q, dim);

      for (int rep = 0; rep < 5; ++rep) {
        dirac.randomiseMatrices(rng);

        const int x = static_cast<int>(rng.getUniformInt(0, dirac.getNumMatrices() - 1));
        const int row = static_cast<int>(rng.getUniformInt(0, dirac.getMatrixDimension() - 1));
        const int col = static_cast<int>(rng.getUniformInt(0, dirac.getMatrixDimension() - 1));
        const arma::cx_double z = randomPerturbation(rng, row, col);

        const double analytic = BarrettGlaserAction::delta2(dirac, x, row, col, z);
        const double direct = calculateDirectDelta2(dirac, x, row, col, z);

        EXPECT_NEAR(direct, analytic, kTolerance)
            << "p=" << p << ", q=" << q << ", dim=" << dim
            << ", matrix=" << x << ", element=(" << row << ", " << col << ")";
      }
    }
  }
}

// 1. Invariant: Action variation under Hermitian rank-1 update Tr(D^4(M + delta_M)) - Tr(D^4(M)).
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Eq. (3.2).
// 3. Dual-Path Verification: BarrettGlaserAction::delta4 vs calculateDirectDelta4.
// 4. Tolerance Rationale: kTolerance (1e-7) bounds double-precision trace cancellation.
TEST(DeltaTests, Delta4MatchesTraceDifference) {
  const StdRng rng(42);

  for (const auto& [p, q] : kGeometries) {
    for (int dim = 2; dim <= 4; ++dim) {
      DiracOperator dirac(p, q, dim);

      for (int rep = 0; rep < 5; ++rep) {
        dirac.randomiseMatrices(rng);

        const int x = static_cast<int>(rng.getUniformInt(0, dirac.getNumMatrices() - 1));
        const int row = static_cast<int>(rng.getUniformInt(0, dirac.getMatrixDimension() - 1));
        const int col = static_cast<int>(rng.getUniformInt(0, dirac.getMatrixDimension() - 1));
        const arma::cx_double z = randomPerturbation(rng, row, col);

        const double analytic = BarrettGlaserAction::delta4(dirac, x, row, col, z);
        const double direct = calculateDirectDelta4(dirac, x, row, col, z);

        EXPECT_NEAR(direct, analytic, kTolerance)
            << "p=" << p << ", q=" << q << ", dim=" << dim
            << ", matrix=" << x << ", element=(" << row << ", " << col << ")";
      }
    }
  }
}

// 1. Invariant: Total action variation Delta S = g_2 * Delta_2 + g_4 * Delta_4.
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Eq. (1.1).
// 3. Dual-Path Verification: BarrettGlaserAction::delta24 vs calculateDirectDelta24.
// 4. Tolerance Rationale: kTolerance (1e-7) bounds double-precision trace cancellation.
TEST(DeltaTests, Delta24MatchesActionDifference) {
  const std::vector<BarrettGlaserAction> actions = {
      BarrettGlaserAction(-3.0, 1.0),
      BarrettGlaserAction(1.5, 0.5)};
  const StdRng rng(42);

  for (const auto& [p, q] : kGeometries) {
    for (int dim = 2; dim <= 4; ++dim) {
      DiracOperator dirac(p, q, dim);

      for (int rep = 0; rep < 5; ++rep) {
        dirac.randomiseMatrices(rng);

        const int x = static_cast<int>(rng.getUniformInt(0, dirac.getNumMatrices() - 1));
        const int row = static_cast<int>(rng.getUniformInt(0, dirac.getMatrixDimension() - 1));
        const int col = static_cast<int>(rng.getUniformInt(0, dirac.getMatrixDimension() - 1));
        const arma::cx_double z = randomPerturbation(rng, row, col);

        for (const auto& action : actions) {
          const double analytic = action.delta24(dirac, x, row, col, z);
          const double direct = calculateDirectDelta24(action, dirac, x, row, col, z);

          EXPECT_NEAR(direct, analytic, kTolerance)
              << "p=" << p << ", q=" << q << ", dim=" << dim
              << ", matrix=" << x << ", element=(" << row << ", " << col << ")"
              << ", g2=" << action.getG2() << ", g4=" << action.getG4();
        }
      }
    }
  }
}

// 1. Invariant: Diagonal variation with purely real z (delta_M = 2 * Re(z) * E_{ii}).
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Section 3.
// 3. Dual-Path Verification: BarrettGlaserAction::delta24 vs calculateDirectDelta24 with fixed seed.
// 4. Tolerance Rationale: kTolerance (1e-7) bounds double-precision trace cancellation.
TEST(DeltaTests, DeterministicDiagonalUpdateMatchesActionDifference) {
  const BarrettGlaserAction action(-2.5, 1.2);
  const StdRng rng(42);
  DiracOperator dirac(1, 1, 4);
  dirac.randomiseMatrices(rng);

  const int x = 0;
  const int row = 1;
  const int col = 1;
  const arma::cx_double z(0.3, 0.0);

  const double analytic = action.delta24(dirac, x, row, col, z);
  const double direct = calculateDirectDelta24(action, dirac, x, row, col, z);

  EXPECT_NEAR(direct, analytic, kTolerance);
}

// 1. Invariant: Off-diagonal variation with complex z (delta_M = z * E_{ij} + conj(z) * E_{ji}).
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Section 3.
// 3. Dual-Path Verification: BarrettGlaserAction::delta24 vs calculateDirectDelta24 with fixed seed.
// 4. Tolerance Rationale: kTolerance (1e-7) bounds double-precision trace cancellation.
TEST(DeltaTests, DeterministicOffDiagonalUpdateMatchesActionDifference) {
  const BarrettGlaserAction action(-2.5, 1.2);
  const StdRng rng(42);
  DiracOperator dirac(1, 1, 4);
  dirac.randomiseMatrices(rng);

  const int x = 0;
  const int row = 1;
  const int col = 2;
  const arma::cx_double z(0.2, -0.4);

  const double analytic = action.delta24(dirac, x, row, col, z);
  const double direct = calculateDirectDelta24(action, dirac, x, row, col, z);

  EXPECT_NEAR(direct, analytic, kTolerance);
}
