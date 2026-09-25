//
// Created for RFL milestone v0.3.0 (Issue #56)
// Public API static contract tests and type trait assurances.
//

#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "IAction.hpp"
#include "IAlgorithm.hpp"
#include "IDiracOperator.hpp"
#include "IRng.hpp"
#include "Simulation.hpp"
#include "StdRng.hpp"
#include <armadillo>
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>
#include <utility>

namespace {

// ============================================================================
// Abstract Interface Invariants
// ============================================================================
static_assert(std::is_abstract_v<IDiracOperator>, "IDiracOperator must be an abstract base interface.");
static_assert(std::has_virtual_destructor_v<IDiracOperator>, "IDiracOperator must define a virtual destructor.");

static_assert(std::is_abstract_v<IAction>, "IAction must be an abstract base interface.");
static_assert(std::has_virtual_destructor_v<IAction>, "IAction must define a virtual destructor.");

static_assert(std::is_abstract_v<IAlgorithm>, "IAlgorithm must be an abstract base interface.");
static_assert(std::has_virtual_destructor_v<IAlgorithm>, "IAlgorithm must define a virtual destructor.");

static_assert(std::is_abstract_v<IRng>, "IRng must be an abstract base interface.");
static_assert(std::has_virtual_destructor_v<IRng>, "IRng must define a virtual destructor.");
static_assert(std::is_same_v<decltype(&IRng::getGaussian), double (IRng::*)(double) const>,
              "IRng::getGaussian must accept sigma, return double, and be const.");
static_assert(std::is_same_v<decltype(&IRng::getUniform), double (IRng::*)() const>,
              "IRng::getUniform must return double and be const.");
static_assert(std::is_same_v<decltype(&IRng::getUniformInt), uint64_t (IRng::*)(uint64_t, uint64_t) const>,
              "IRng::getUniformInt must accept (min, max), return uint64_t, and be const.");

// ============================================================================
// DiracOperator Contract Assertions
// ============================================================================
static_assert(std::is_base_of_v<IDiracOperator, DiracOperator>, "DiracOperator must inherit from IDiracOperator.");
static_assert(std::is_final_v<DiracOperator>, "DiracOperator must be declared final.");

// Constructibility and value semantics
static_assert(!std::is_default_constructible_v<DiracOperator>, "DiracOperator default constructor must not be available.");
static_assert(std::is_constructible_v<DiracOperator, int, int, int>, "DiracOperator must be constructible from (p, q, dim).");
static_assert(std::is_copy_constructible_v<DiracOperator>, "DiracOperator must be copy-constructible.");
static_assert(std::is_copy_assignable_v<DiracOperator>, "DiracOperator must be copy-assignable.");
static_assert(std::is_move_constructible_v<DiracOperator>, "DiracOperator must be move-constructible.");
static_assert(std::is_move_assignable_v<DiracOperator>, "DiracOperator must be move-assignable.");

// Core public inspector signatures and const-qualifications
static_assert(std::is_same_v<decltype(&DiracOperator::getType), std::pair<int, int> (DiracOperator::*)() const>,
              "DiracOperator::getType must return std::pair<int, int> and be const.");

static_assert(std::is_same_v<decltype(&DiracOperator::getMatrixDimension), int (DiracOperator::*)() const>,
              "DiracOperator::getMatrixDimension must return int and be const.");

static_assert(std::is_same_v<decltype(&DiracOperator::getEigenvalues), arma::vec (DiracOperator::*)() const>,
              "DiracOperator::getEigenvalues must return arma::vec and be const.");

static_assert(std::is_same_v<decltype(&DiracOperator::getDiracMatrix), arma::cx_mat (DiracOperator::*)() const>,
              "DiracOperator::getDiracMatrix must return arma::cx_mat and be const.");

static_assert(std::is_same_v<decltype(&DiracOperator::randomiseMatrices), void (DiracOperator::*)(const IRng&)>,
              "DiracOperator::randomiseMatrices must be non-const.");

using ConstGetMatrices = const std::vector<arma::cx_mat>& (DiracOperator::*)() const;
using MutGetMatrices = std::vector<arma::cx_mat>& (DiracOperator::*)();
static_assert(std::is_same_v<ConstGetMatrices, decltype(static_cast<ConstGetMatrices>(&DiracOperator::getMatrices))>,
              "DiracOperator must provide a const getMatrices() overload.");
static_assert(std::is_same_v<MutGetMatrices, decltype(static_cast<MutGetMatrices>(&DiracOperator::getMatrices))>,
              "DiracOperator must provide a mutable getMatrices() overload.");

// ============================================================================
// Action Contract Assertions
// ============================================================================
static_assert(std::is_base_of_v<IAction, Action>, "Action must inherit from IAction.");
static_assert(std::is_default_constructible_v<Action>, "Action must be default-constructible.");
static_assert(std::is_constructible_v<Action, double, double>, "Action must be constructible from (g_2, g_4).");
static_assert(std::is_constructible_v<Action, double>, "Action must be constructible from single coupling (g_2).");

// Value type semantics
static_assert(std::is_copy_constructible_v<Action>, "Action must be copy-constructible.");
static_assert(std::is_copy_assignable_v<Action>, "Action must be copy-assignable.");
static_assert(std::is_move_constructible_v<Action>, "Action must be move-constructible.");
static_assert(std::is_move_assignable_v<Action>, "Action must be move-assignable.");

// Method signatures and const-qualifications
static_assert(std::is_same_v<decltype(&Action::getG2), double (Action::*)() const>,
              "Action::getG2 must return double and be const.");

static_assert(std::is_same_v<decltype(&Action::getG4), double (Action::*)() const>,
              "Action::getG4 must return double and be const.");

static_assert(std::is_same_v<decltype(&Action::setParams), void (Action::*)(double, double)>,
              "Action::setParams must accept (double, double) and return void.");

static_assert(std::is_same_v<decltype(&Action::calculateS), double (Action::*)(const IDiracOperator&) const>,
              "Action::calculateS must accept const IDiracOperator&, return double, and be const.");

// ============================================================================
// StdRng Contract Assertions
// ============================================================================
static_assert(std::is_base_of_v<IRng, StdRng>, "StdRng must inherit from IRng.");
static_assert(std::is_final_v<StdRng>, "StdRng must be declared final.");
static_assert(std::is_default_constructible_v<StdRng>, "StdRng must be default-constructible.");
static_assert(std::is_constructible_v<StdRng, const uint64_t>, "StdRng must be constructible from 64-bit seed.");

static_assert(std::is_same_v<decltype(&StdRng::getGaussian), double (StdRng::*)(const double) const>,
              "StdRng::getGaussian must accept sigma, return double, and be const.");

static_assert(std::is_same_v<decltype(&StdRng::getUniform), double (StdRng::*)() const>,
              "StdRng::getUniform must return double and be const.");

static_assert(std::is_same_v<decltype(&StdRng::getUniformInt), uint64_t (StdRng::*)(const uint64_t, const uint64_t) const>,
              "StdRng::getUniformInt must accept (min, max), return uint64_t, and be const.");

// ============================================================================
// Metropolis Contract Assertions
// ============================================================================
static_assert(std::is_base_of_v<IAlgorithm, Metropolis>, "Metropolis must inherit from IAlgorithm.");
static_assert(std::is_final_v<Metropolis>, "Metropolis must be declared final.");
static_assert(!std::is_default_constructible_v<Metropolis>, "Metropolis default constructor must be deleted.");
static_assert(std::is_constructible_v<Metropolis, std::unique_ptr<Action>&&, const double, const int, std::unique_ptr<IRng>&&>,
              "Metropolis must be constructible from (action, scale, steps, rng).");
static_assert(!std::is_copy_constructible_v<Metropolis>, "Metropolis must be non-copyable.");
static_assert(!std::is_copy_assignable_v<Metropolis>, "Metropolis must not be copy-assignable.");

static_assert(std::is_same_v<decltype(&Metropolis::updateDirac), double (Metropolis::*)(IDiracOperator&) const>,
              "Metropolis::updateDirac must accept IDiracOperator&, return double, and be const.");

// ============================================================================
// Simulation Contract Assertions
// ============================================================================
static_assert(!std::is_default_constructible_v<Simulation>, "Simulation default constructor must be deleted.");
static_assert(std::is_constructible_v<Simulation, std::unique_ptr<DiracOperator>&&, std::unique_ptr<IAlgorithm>&&>,
              "Simulation must be constructible from (dirac, algorithm).");
static_assert(!std::is_copy_constructible_v<Simulation>, "Simulation must be non-copyable.");
static_assert(!std::is_copy_assignable_v<Simulation>, "Simulation must not be copy-assignable.");

static_assert(std::is_same_v<decltype(&Simulation::run), double (Simulation::*)() const>,
              "Simulation::run must return double and be const.");

static_assert(std::is_same_v<decltype(&Simulation::getDiracOperator), const DiracOperator& (Simulation::*)() const>,
              "Simulation::getDiracOperator must return const DiracOperator& and be const.");

}// namespace

// ============================================================================
// Runtime Sanity Tests
// ============================================================================
TEST(ApiContractTest, EnduringPublicContractSanity) {
  // Verify DiracOperator runtime instantiation and inspector contracts.
  DiracOperator dirac(1, 3, 6);
  EXPECT_EQ(dirac.getType().first, 1);
  EXPECT_EQ(dirac.getType().second, 3);
  EXPECT_EQ(dirac.getMatrixDimension(), 6);

  arma::vec evs = dirac.getEigenvalues();
  EXPECT_GT(evs.n_elem, 0U);

  // Verify Action runtime instantiation and contract.
  Action action(-2.0, 1.0);
  EXPECT_DOUBLE_EQ(action.getG2(), -2.0);
  EXPECT_DOUBLE_EQ(action.getG4(), 1.0);

  double s = action.calculateS(dirac);
  EXPECT_TRUE(std::isfinite(s));

  // Verify StdRng runtime instantiation.
  StdRng rng(42UL);
  double u = rng.getUniform();
  EXPECT_GE(u, 0.0);
  EXPECT_LT(u, 1.0);

  // Verify Metropolis and Simulation assembly.
  auto act_ptr = std::make_unique<Action>(-2.0, 1.0);
  auto rng_ptr = std::make_unique<StdRng>(42UL);
  auto metro = std::make_unique<Metropolis>(std::move(act_ptr), 0.1, 1, std::move(rng_ptr));
  auto dirac_ptr = std::make_unique<DiracOperator>(1, 3, 6);

  Simulation sim(std::move(dirac_ptr), std::move(metro));
  EXPECT_EQ(sim.getDiracOperator().getMatrixDimension(), 6);

  double rate = sim.run();
  EXPECT_GE(rate, 0.0);
  EXPECT_LE(rate, 1.0);
}
