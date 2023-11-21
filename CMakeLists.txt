cmake_minimum_required(VERSION 3.7...3.23)

if(${CMAKE_VERSION} VERSION_LESS 3.23)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.23)
endif()

project(
        RFL
        VERSION 0.5
        DESCRIPTION
        "Random Fuzzy Library (RFL) - A library for Markov Chain Monte Carlo simulations of Finite Random Non-commutative geometries."
        LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  # Nicely supports folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  enable_testing()

  add_compile_options(-Wall -Wextra -Wpedantic)
  add_compile_options(-fsanitize=address )
  add_link_options(-fsanitize=address )
endif()

# Include ccache to speed up rebuild times.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Ccache)

option(${PROJECT_NAME}_ENABLE_CODE_COVERAGE "Enable code coverage through GCC." ON)


include(Doxygen)
doxygen_add_docs(doxygen
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/RFL_source/new_source
        ${PROJECT_SOURCE_DIR}/examples
        COMMENT "Generate documentation for entire project")

# GSL is included with Armadillo above, but we need more.
find_package(Armadillo REQUIRED)
find_package(GSL REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
add_subdirectory(vendors)
include(GoogleTest)

add_subdirectory(RFL_source)
add_subdirectory(examples)

# # Generate coverage using gcov
# if (${CMAKE_PROJECT_NAME}_ENABLE_CODE_COVERAGE)
#   target_compile_options(${CMAKE_PROJECT_NAME} PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
#   target_link_options(${CMAKE_PROJECT_NAME} PUBLIC -fprofile-arcs -ftest-coverage)
# endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


# AN DEMO/EXPERIMENT AREA
add_subdirectory(playground)
