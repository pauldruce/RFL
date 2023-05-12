add_executable(
        new_RFLTests
        tAction.cpp tClifford.cpp tDiracOperator.cpp tHamiltonian.cpp tMetropolis.cpp tGslRng.cpp tSimulation.cpp)

target_include_directories(new_RFLTests PUBLIC ${new_RFL_INCLUDE_DIR})
target_link_libraries(new_RFLTests gtest_main gmock_main new_RFL ${ARMADILLO_LIBRARIES})

gtest_discover_tests(new_RFLTests DISCOVERY_TIMEOUT 600)


add_executable(BenchmarkTests
        tBenchmark.cpp)
target_include_directories(BenchmarkTests PUBLIC ${new_RFL_INCLUDE_DIR})
target_include_directories(BenchmarkTests PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(BenchmarkTests gtest_main RFL new_RFL ${ARMADILLO_LIBRARIES})

gtest_discover_tests(BenchmarkTests DISCOVERY_TIMEOUT 600)
