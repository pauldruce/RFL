set(INCLUDES ${INCLUDE_DIR}/clifford.hpp ${INCLUDE_DIR}/geometry.hpp)

add_library(RFL)

target_sources(
  RFL
  PRIVATE action.cpp
          clifford.cpp
          derivative.cpp
          hamiltonian.cpp
          HMC.cpp
          HMC_core.cpp
          leapfrog.cpp
          metropolis.cpp
          misc.cpp
          MMC.cpp
          MMC_core.cpp
  PUBLIC ${INCLUDES})

target_include_directories(RFL PUBLIC ${INCLUDE_DIR})
target_link_libraries(RFL ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})