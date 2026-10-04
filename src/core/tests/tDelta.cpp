//
// Unit tests verifying analytic trace variations against direct action differences.
//

#include "../BarrettGlaser/BarrettGlaserAction.hpp"
#include "../DiracOperator.hpp"
#include "../StdRng.hpp"
#include <algorithm>
#include <armadillo>
#include <cmath>
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
    for (int dim = 1; dim <= 4; ++dim) {
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
    for (int dim = 1; dim <= 4; ++dim) {
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
    for (int dim = 1; dim <= 4; ++dim) {
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

// 1. Invariant: Action derivative along Hermitian tangent direction matches directional gradient.
//    lim_{eps -> 0} (S(M_k + eps * H) - S(M_k - eps * H)) / (2 * eps) == Re(Tr(nabla_{M_k} S * H)).
// 2. Reference: EP-4 Computational Physics Verification Suite, Section 3.
// 3. Dual-Path Verification: Analytical derDirac24 vs central finite difference of calculateS.
// 4. Tolerance Rationale: Central stencil step eps = 1e-6 balances truncation error O(eps^2) and
//    floating-point cancellation O(eps_mach / eps) to optimal O(eps_mach^(2/3)) precision.
TEST(DeltaTests, GradientMatchesFiniteDifference) {
  const std::vector<BarrettGlaserAction> actions = {
      BarrettGlaserAction(-2.0, 1.5),
      BarrettGlaserAction(1.0, 0.5)};
  const StdRng rng(42);
  constexpr double kEps = 1e-6;
  constexpr double kGradientRelTol = 1e-7;

  const std::vector<std::pair<int, int>> test_geometries = {
      {1, 0},
      {0, 1},
      {1, 1},
      {2, 0},
      {2, 2}};

  for (const auto& [p, q] : test_geometries) {
    for (int dim = 2; dim <= 3; ++dim) {
      DiracOperator dirac(p, q, dim);
      dirac.randomiseMatrices(rng);

      for (const auto& action : actions) {
        for (int k = 0; k < dirac.getNumMatrices(); ++k) {
          // Construct random Hermitian perturbation direction H = 0.5 * (W + W†).
          arma::cx_mat direction(dim, dim);
          for (int r = 0; r < dim; ++r) {
            for (int c = 0; c < dim; ++c) {
              direction(r, c) = arma::cx_double(-0.5 + rng.getUniform(), -0.5 + rng.getUniform());
            }
          }
          direction = 0.5 * (direction + direction.t());

          DiracOperator dirac_plus = dirac;
          dirac_plus.getMatrix(k) += kEps * direction;

          DiracOperator dirac_minus = dirac;
          dirac_minus.getMatrix(k) -= kEps * direction;

          const double s_plus = action.calculateS(dirac_plus);
          const double s_minus = action.calculateS(dirac_minus);
          const double numerical_derivative = (s_plus - s_minus) / (2.0 * kEps);

          const arma::cx_mat grad_herm = action.derDirac24(dirac, k, true);
          const arma::cx_mat grad_raw = action.derDirac24(dirac, k, false);
          const double analytical_derivative = std::real(arma::trace(grad_herm * direction));
          const double raw_derivative = std::real(arma::trace(grad_raw * direction));

          const double scale =
              std::max({1.0, std::abs(analytical_derivative), 1e-4 * (std::abs(s_plus) + std::abs(s_minus))});

          // For any Hermitian perturbation direction, projected and raw variations match:
          EXPECT_NEAR(raw_derivative, analytical_derivative, 1e-10 * scale);

          EXPECT_NEAR(numerical_derivative, analytical_derivative, kGradientRelTol * scale)
              << "p=" << p << ", q=" << q << ", dim=" << dim << ", matrix=" << k
              << ", g2=" << action.getG2() << ", g4=" << action.getG4();
        }
      }
    }
  }
}

// 1. Invariant: Linear scaling of action variation and gradients with coupling constants.
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Eq. (1.1).
// 3. Invariant: delta24(0, g_4) == g_4 * delta4, and derDirac24(0, g_4) == g_4 * derDirac4.
// 4. Tolerance Rationale: Exact double-precision floating-point linearity within machine precision.
TEST(DeltaTests, QuarticCouplingScalesLinearly) {
  const StdRng rng(42);
  DiracOperator dirac(1, 1, 3);
  dirac.randomiseMatrices(rng);

  const double delta_2 = 1.7;
  const double delta_4 = -2.3;

  const BarrettGlaserAction action_unit(0.0, 1.0);
  const BarrettGlaserAction action_scaled(0.0, 2.5);

  EXPECT_DOUBLE_EQ(action_unit.delta24(delta_2, delta_4), delta_4);
  EXPECT_DOUBLE_EQ(action_scaled.delta24(delta_2, delta_4), 2.5 * delta_4);

  const arma::cx_mat grad_unit = action_unit.derDirac24(dirac, 0, true);
  const arma::cx_mat grad_scaled = action_scaled.derDirac24(dirac, 0, true);

  const double max_diff = arma::max(arma::max(arma::abs(grad_scaled - 2.5 * grad_unit)));
  EXPECT_NEAR(max_diff, 0.0, 1e-12);
}

// 1. Invariant: Action invariance under Dirac field rescaling (EP-4, REQ-011).
//    S(alpha * D; alpha^(-2) * g_2, alpha^(-4) * g_4) == S(D; g_2, g_4).
//    Specifically for alpha = g_4^(1/4): S(g_4^(1/4) * D; g_2 / sqrt(g_4), 1.0) == S(D; g_2, g_4).
// 2. Reference: Barrett & Glaser (2016), arXiv:1510.01377, Section 2; docs/theory/barrett_glaser.md#field-rescaling--literature-conventions.
// 3. Tolerance Rationale: Exact mathematical identity bounded by floating-point arithmetic precision (1e-11).
TEST(DeltaTests, FieldRescalingEquivalence) {
  const StdRng rng(42);
  const double g_2 = -2.5;
  const double g_4 = 4.0;
  const BarrettGlaserAction physical_action(g_2, g_4);

  // Field rescaling factor alpha = g_4^(1/4) = 4^(1/4) = sqrt(2).
  const double alpha = std::pow(g_4, 0.25);
  const double g_2_tilde = g_2 / std::sqrt(g_4);
  const double g_4_tilde = 1.0;
  const BarrettGlaserAction rescaled_action(g_2_tilde, g_4_tilde);

  for (const auto& [p, q] : kGeometries) {
    DiracOperator dirac(p, q, 3);
    dirac.randomiseMatrices(rng);

    // Create rescaled Dirac operator D_tilde = alpha * D
    DiracOperator dirac_rescaled = dirac;
    for (int k = 0; k < dirac.getNumMatrices(); ++k) {
      dirac_rescaled.getMatrix(k) *= alpha;
    }

    const double s_physical = physical_action.calculateS(dirac);
    const double s_rescaled = rescaled_action.calculateS(dirac_rescaled);

    EXPECT_NEAR(s_physical, s_rescaled, 1e-11 * std::max(1.0, std::abs(s_physical)))
        << "p=" << p << ", q=" << q << ", physical=" << s_physical << ", rescaled=" << s_rescaled;
  }
}
