- [ ] Documentation
  - [x] Add readmes everywhere to guidance users through source code, tests and examples.
  - [ ] Add documentation string for each class and public method
  - [ ] Add documentation generation to github actions
- [ ] Transfer methods that can be made static in Action to more appropriate place
- [ ] Investigate multithreading to improve performance.
  - It might be possible to use Armadillo here, rather than do the multithreading myself.
- [ ] Replace raw pointer use to use smart pointers and references.
  - [x] DiracOperator
  - [x] Action
  - [ ] Metropolis
  - [ ] Hamiltonian
- [ ] Refactor the github actions to have reusable setup actions for linux, windows, macos.
- [ ] Create github actions to create built libraries for macOS, Windows and linux OS.
- [ ] Make an interface for DiracOperator, rather than using concrete DiracOperator in
  Simulation class and Algorithm.

## Completed

- [x] Unit tests for:
  - [x] Hamiltonian Class
  - [x] Metropolis Class
  - [x] Simulation
- [x] Create an Action interface/abstract class.
- [x] Add a clang-format step to GitHub actions
- [x] Create some demo apps to showcase how to create and run a simple simulation with the new library
- [x] Create an interface for the random number generator and wrap the gsl methods we use in a class.
  - [x] Interface made
  - [x] GSL class created
  - [x] Refactor code to make use of interface rather than GSL directly