add_executable(
        new_RFLTests
        tAction.cpp tClifford.cpp tDiracOperator.cpp tHamiltonian.cpp tMetropolis.cpp tGslRng.cpp tSimulation.cpp)

target_include_directories(new_RFLTests PUBLIC ${new_RFL_INCLUDE_DIR})
target_link_libraries(new_RFLTests gtest_main gmock_main new_RFL ${ARMADILLO_LIBRARIES})

gtest_discover_tests(new_RFLTests DISCOVERY_TIMEOUT 600)
add_subdirectory(performance)