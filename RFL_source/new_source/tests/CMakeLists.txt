add_executable(
        new_RFLTests
        tAction.cpp
        tClifford.cpp tDiracOperator.cpp)

target_include_directories(new_RFLTests PUBLIC ${new_RFL_INCLUDE_DIR})
target_link_libraries(new_RFLTests gtest_main new_RFL ${ARMADILLO_LIBRARIES})
target_compile_options(new_RFLTests PRIVATE -fsanitize=address)
target_link_options(new_RFLTests PRIVATE -fsanitize=address)

include(GoogleTest)
gtest_discover_tests(new_RFLTests)
