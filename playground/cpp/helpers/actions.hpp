// FILE: actions.hpp

#ifndef ACTIONS_HPP
#define ACTIONS_HPP

#include <string>

// Abstract base class representing the Strategy interface
class Action {
public:
  virtual ~Action() = default;
  // Pure virtual function to be implemented by concrete strategies
  virtual std::string calculate() = 0;
};

// Concrete Strategy implementing the Action interface
class quadraticAction final : public Action {
  std::string calculate() override {
    return "calculated by quadratic action";
  }
};

// Another Concrete Strategy implementing the Action interface
class barrettGlaserAction final : public Action {
  std::string calculate() override {
    return "calculated by Barrett-Glaser action";
  }
};

#endif// ACTIONS_HPP