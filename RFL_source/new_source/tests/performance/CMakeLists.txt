add_executable(new_RFL_PerformanceTests
        tBenchmark.cpp)
target_include_directories(new_RFL_PerformanceTests PUBLIC ${new_RFL_INCLUDE_DIR})
target_include_directories(new_RFL_PerformanceTests PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(new_RFL_PerformanceTests gtest_main RFL new_RFL ${ARMADILLO_LIBRARIES})

gtest_discover_tests(new_RFL_PerformanceTests DISCOVERY_TIMEOUT 600)
