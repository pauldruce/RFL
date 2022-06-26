add_executable(RFLTests_doctest
        doctest_main.cpp
        action_tests.cpp
        compilation_test.cpp
        delta_tests.cpp)

target_include_directories(RFLTests_doctest PUBLIC ${DOCTEST_INCLUDE_DIR}
                                           ${INCLUDE_DIR})
target_link_libraries(RFLTests_doctest RFL ${ARMADILLO_LIBRARIES})
target_compile_options(RFLTests_doctest PRIVATE -fsanitize=address)
target_link_options(RFLTests_doctest PRIVATE -fsanitize=address)

# Use CMake module to discover Doctest tests.
include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(RFLTests_doctest)
