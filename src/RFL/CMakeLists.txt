cmake_minimum_required(VERSION 3.20...3.31)


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
        if(NOT SKBUILD)
            add_compile_options(-fsanitize=address)
            add_link_options(-fsanitize=address)
        endif()
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
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND NOT SKBUILD)
        add_subdirectory(${RFL_SOURCE_DIR}/tests)
endif()

# Python Bindings
option(BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" ON)
if(BUILD_PYTHON_BINDINGS)
    include(FetchContent)

    # Fetch carma (NumPy/Armadillo bridge)
    FetchContent_Declare(
        carma
        GIT_REPOSITORY https://github.com/RUrlus/carma.git
        GIT_TAG        v0.8.0
    )
    FetchContent_Populate(carma)
    add_library(carma::carma INTERFACE IMPORTED)
    set_target_properties(carma::carma PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${carma_SOURCE_DIR}/include"
    )

    # Fetch pybind11
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.1
    )
    FetchContent_MakeAvailable(pybind11)

    add_subdirectory(python_bindings)
endif()