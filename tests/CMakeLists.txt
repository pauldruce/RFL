add_executable(RFLTests    doctest_main.cpp
                            action_tests.cpp)

target_include_directories(
        RFLTests
        PUBLIC  ${DOCTEST_INCLUDE_DIR}
                ${INCLUDE_DIR})
target_link_libraries(RFLTests RFL)