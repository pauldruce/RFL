execute_process(
    COMMAND git rev-parse --show-toplevel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_PROJECT_ROOT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)