//
// Created by Paul Druce on 26/03/2023.
//

#include <iostream>
#include <memory>

class Action {
public:
  virtual ~Action() = default;
  virtual std::string calculate() = 0;
};

class quadraticAction final : public Action {
  std::string calculate() override {
    return "calculated by quadratic action";
  };
};

class barrettGlaserAction final : public Action {
  std::string calculate() override {
    return "calculated by Barrett-Glaser action";
  }
};

class ActionManager {
public:
  explicit ActionManager(std::unique_ptr<Action>&& input_action = {})
      : m_action(std::move(input_action)) {
  }

  void SetAction(std::unique_ptr<Action>&& new_action) {
    m_action = std::move(new_action);
  }

  void Print() const {
    const auto action_val = m_action->calculate();
    std::cout << action_val << std::endl;
  }

private:
  std::unique_ptr<Action> m_action;
};

int main() {
  ActionManager AM(std::make_unique<quadraticAction>());
  AM.Print();

  AM.SetAction(std::make_unique<barrettGlaserAction>());
  AM.Print();
}