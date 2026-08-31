add_executable(
        rfl_legacy_tests action_tests.cpp compilation_tests.cpp delta_tests.cpp
        derivative_tests.cpp geometry.test.cpp clifford.test.cpp)

target_include_directories(rfl_legacy_tests PUBLIC ${legacy_RFL_INCLUDE_DIR})
target_link_libraries(rfl_legacy_tests PRIVATE gtest_main rfl_legacy ${ARMADILLO_LIBRARIES})

gtest_discover_tests(rfl_legacy_tests DISCOVERY_TIMEOUT 600)
