//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_ACTION_CPP_H
#define RFL_ACTION_CPP_H
#include <armadillo>
#include <string>
#include <gsl/gsl_rng.h>
#include "DiracOperator.hpp"

class Action {
 public:
  // CONSTRUCTORS AND DESTRUCTORS
  Action(double g2, double g4);
  Action() : g2(0.0), g4(0.0) {};
  ~Action() = default;

  // METHODS
 public:
  void set_g2(double value);
  void set_g4(double value);
  void set_params(double g2, double g4);
  double get_g2() const { return g2; }
  double get_g4() const { return g4; }

 public:
  double calculate_S(const DiracOperator &D) const;
  double calculate_S_from_dirac(const DiracOperator &D) const;
  double dirac2(const DiracOperator &D) const;
  double dirac4(const DiracOperator &D) const;

  void print_S(const DiracOperator &D, std::ostream &out) const {
	out << dirac2(D) << " " << dirac4(D) << std::endl;
  }

 private:
  double g2, g4;

  double compute_A4(const DiracOperator &D,
					const int &i1,
					const int &i2,
					const int &i3,
					const int &i4) const;
  double compute_A2(const DiracOperator &D, const int &i1, const int &i2) const;
  double compute_A(const DiracOperator &D, const int &i) const;
};
#endif // RFL_ACTION_CPP_H
