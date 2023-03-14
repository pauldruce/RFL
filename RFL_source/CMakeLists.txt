# Build RFL library
set(RFL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(RFL_INCLUDE_DIR ${RFL_SOURCE_DIR}/include)
add_subdirectory(${RFL_SOURCE_DIR}/source)

# Add new version of library
set(new_RFL_SOURCE_DIR "${RFL_SOURCE_DIR}/new_source")
set(new_RFL_INCLUDE_DIR "${RFL_SOURCE_DIR}/new_source")
add_subdirectory(${RFL_SOURCE_DIR}/new_source)

# Add tests if building tests.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(${RFL_SOURCE_DIR}/tests)
endif ()