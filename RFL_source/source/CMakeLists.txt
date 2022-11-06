set(RFL_INCLUDES ${RFL_INCLUDE_DIR}/Cliff.hpp ${RFL_INCLUDE_DIR}/Geom24.hpp)

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
  PUBLIC ${RFL_INCLUDES})

target_include_directories(RFL PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(RFL ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})
