add_executable(
        RFLTests_GTest action_tests.cpp compilation_tests.cpp delta_tests.cpp
        derivative_tests.cpp geometry.test.cpp clifford.test.cpp)

target_include_directories(RFLTests_GTest PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(RFLTests_GTest gtest_main RFL ${ARMADILLO_LIBRARIES})

gtest_discover_tests(RFLTests_GTest DISCOVERY_TIMEOUT 600)
