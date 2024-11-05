#include "helpers/actions.hpp"
#include <iostream>
#include <memory>
#include <string>

// Context class that uses an Action strategy
class ActionManager {
public:
  explicit ActionManager(std::unique_ptr<Action>&& input_action = {})
      : m_action(std::move(input_action)) {
  }

  // Method to set a new Action strategy
  void SetAction(std::unique_ptr<Action>&& new_action) {
    m_action = std::move(new_action);
  }

  // Method to print the result of the current Action's calculate method
  void DisplayActionValue() const {
    std::cout << "Calculate the action value...\n";
    const auto action_val = m_action->calculate();
    std::cout << "The action value is: \"" << action_val << "\"" << std::endl;
  }

private:
  std::unique_ptr<Action> m_action;// Pointer to the current Action strategy
};

int main() {
  std::cout << "This example demonstrates the Strategy Design Pattern.\n";

  std::cout << "Creating ActionManager with quadraticAction...\n";
  ActionManager actionManager(std::make_unique<quadraticAction>());
  actionManager.DisplayActionValue();// Output: calculated by quadratic action

  std::cout << "\nChanging Action to barrettGlaserAction...\n";
  actionManager.SetAction(std::make_unique<barrettGlaserAction>());
  actionManager.DisplayActionValue();// Output: calculated by Barrett-Glaser action

  std::cout << "\nThe ActionManager can change its behavior at runtime by switching between different Action strategies.";

  return 0;
}