- [ ] Unit tests for:
  - [x] Hamiltonian Class
  - [ ] Metropolis Class
  - [ ] Simulation -> likely an end-to-end test of some description.
    - Look into how we should test the random-ness. How do we set a seed for this in a way that we can use in testing.
- [ ] Create an Action interface/abstract class. 

- [ ] Transfer methods that can be made static in Action to more appropriate place
- [ ] Add a clang-format step to GitHub actions
  - [ ] Can we add it as a step in CMake build?
- [ ] Investigate the use of multithreading to improve performance. 
  - Potentially can just use Armadillo here, rather than do the multithreading myself.
- [ ] Modify the benchmark examples to average over a few runs. 


## Completed
- [x] Create some demo apps to showcase how to create and run a simple simulation with the new library
- [x] Create an interface for the random number generator and wrap the gsl methods we use in a class.
  - [x] Interface made
  - [x] GSL class created
  - [x] Refactor code to make use of interface rather than GSL directly