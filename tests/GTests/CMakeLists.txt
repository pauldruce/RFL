add_executable(RFLTests_GTest
        action_tests.cpp compilation_tests.cpp delta_tests.cpp derivative_tests.cpp)

target_include_directories(RFLTests_GTest PUBLIC
        ${INCLUDE_DIR})
target_link_libraries(RFLTests_GTest gtest_main RFL ${ARMADILLO_LIBRARIES})
target_compile_options(RFLTests_GTest PRIVATE -fsanitize=address)
target_link_options(RFLTests_GTest PRIVATE -fsanitize=address)

include(GoogleTest)
gtest_discover_tests(RFLTests_GTest)