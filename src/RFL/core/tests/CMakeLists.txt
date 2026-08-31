add_executable(
        rfl_tests
        tAction.cpp tClifford.cpp tDiracOperator.cpp tHamiltonian.cpp tMetropolis.cpp tGslRng.cpp tSimulation.cpp)

target_include_directories(rfl_tests PUBLIC ${core_RFL_INCLUDE_DIR})
target_link_libraries(rfl_tests PRIVATE gtest_main gmock_main rfl_core ${ARMADILLO_LIBRARIES})

gtest_discover_tests(rfl_tests DISCOVERY_TIMEOUT 600)
add_subdirectory(performance)