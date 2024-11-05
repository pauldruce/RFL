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


- [ ] Change new_RFL source code to use an include folder pattern, rather than
  mix and match the srcs and headers. It's just extra work for no gain.