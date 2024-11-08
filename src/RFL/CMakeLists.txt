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
        add_compile_options(-fsanitize=address)
        add_link_options(-fsanitize=address)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Include ccache to speed up rebuild times.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Ccache)

option(${PROJECT_NAME}_ENABLE_CODE_COVERAGE "Enable code coverage through GCC." ON)

# Include Please's output directory
include(GitRoot)

include(Doxygen)
doxygen_add_docs(doxygen
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/new_RFL
        ${PROJECT_SOURCE_DIR}/examples
        COMMENT "Generate documentation for entire project")

# GSL is included with Armadillo above, but we need more.
include(Armadillo)
include(GSL)
include(LoadGoogleTest)

# Build RFL library
set(RFL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/old_RFL")
set(RFL_INCLUDE_DIR "${RFL_SOURCE_DIR}/include")
add_subdirectory(${RFL_SOURCE_DIR}/source)

# Add new version of library
set(new_RFL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/new_RFL")
set(new_RFL_INCLUDE_DIR "${new_RFL_SOURCE_DIR}")
add_subdirectory(${new_RFL_SOURCE_DIR})

# Add tests if building tests.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
        add_subdirectory(${RFL_SOURCE_DIR}/tests)
endif()