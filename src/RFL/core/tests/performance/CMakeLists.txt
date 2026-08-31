add_executable(rfl_performance_tests
        tBenchmark.cpp)
target_include_directories(rfl_performance_tests PUBLIC ${core_RFL_INCLUDE_DIR})
target_include_directories(rfl_performance_tests PUBLIC ${legacy_RFL_INCLUDE_DIR})
target_link_libraries(rfl_performance_tests PRIVATE gtest_main rfl_legacy rfl_core ${ARMADILLO_LIBRARIES})

gtest_discover_tests(rfl_performance_tests DISCOVERY_TIMEOUT 600)
