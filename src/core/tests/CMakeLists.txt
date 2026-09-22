add_executable(
        rfl_tests
        tAction.cpp tApiContract.cpp tClifford.cpp tDiracOperator.cpp tHamiltonian.cpp tMetropolis.cpp tStdRng.cpp tSimulation.cpp)

target_include_directories(rfl_tests PUBLIC ${core_RFL_INCLUDE_DIR})
target_link_libraries(rfl_tests PRIVATE gtest_main gmock_main rfl_core ${ARMADILLO_LIBRARIES})

gtest_discover_tests(rfl_tests DISCOVERY_TIMEOUT 600)

if(TARGET rfl_legacy)
    add_subdirectory(performance)
else()
    message(STATUS "rfl_legacy not available: skipping rfl_performance_tests benchmark")
endif()
