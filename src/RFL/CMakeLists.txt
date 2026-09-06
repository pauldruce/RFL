cmake_minimum_required(VERSION 3.20...3.31)


if(DEFINED SKBUILD_PROJECT_VERSION)
    set(RFL_PROJECT_VERSION ${SKBUILD_PROJECT_VERSION})
else()
    set(RFL_PROJECT_VERSION 0.1.0)
endif()

# Auto-detect vcpkg toolchain file on Windows if not explicitly specified.
if(WIN32 AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(_VCPKG_ROOT "")
    if(DEFINED ENV{VCPKG_ROOT})
        set(_VCPKG_ROOT "$ENV{VCPKG_ROOT}")
    elseif(DEFINED ENV{VCPKG_INSTALLATION_ROOT})
        set(_VCPKG_ROOT "$ENV{VCPKG_INSTALLATION_ROOT}")
    endif()

    set(_VCPKG_TOOLCHAIN "${_VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
    if(EXISTS "${_VCPKG_TOOLCHAIN}")
        set(CMAKE_TOOLCHAIN_FILE "${_VCPKG_TOOLCHAIN}" CACHE STRING "vcpkg toolchain file")
    endif()
    unset(_VCPKG_ROOT)
    unset(_VCPKG_TOOLCHAIN)
endif()

project(
        RFL
        VERSION ${RFL_PROJECT_VERSION}
        DESCRIPTION
        "Random Fuzzy Library (RFL) - A library for Markov Chain Monte Carlo simulations of Finite Random Non-commutative geometries."
        LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
        set(CMAKE_CXX_STANDARD 17)
        set(CMAKE_CXX_STANDARD_REQUIRED ON)
        set(CMAKE_CXX_EXTENSIONS OFF)

        # Nicely supports folders in IDEs
        set_property(GLOBAL PROPERTY USE_FOLDERS ON)
        enable_testing()

        if(MSVC)
            add_compile_options(/W4 /utf-8 /permissive-)
            add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN)
        else()
            add_compile_options(-Wall -Wextra -Wpedantic)
            if(NOT SKBUILD)
                add_compile_options(-fsanitize=address)
                add_link_options(-fsanitize=address)
            endif()
        endif()
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Include ccache to speed up rebuild times.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Ccache)

option(${PROJECT_NAME}_ENABLE_CODE_COVERAGE "Enable code coverage through GCC." ON)
option(RFL_FETCH_ARMADILLO "Fetch Armadillo via CMake FetchContent" OFF)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND NOT SKBUILD)
    include(Doxygen)
    if(DOXYGEN_FOUND)
        doxygen_add_docs(doxygen
                ${PROJECT_SOURCE_DIR}/README.md
                ${PROJECT_SOURCE_DIR}/core
                ${PROJECT_SOURCE_DIR}/examples
                COMMENT "Generate documentation for entire project")
    endif()
endif()

# GSL is included with Armadillo above, but we need more.
include(Armadillo)
include(GSL)
if(NOT SKBUILD)
    include(LoadGoogleTest)
endif()

# Build legacy RFL library
set(legacy_RFL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/legacy")
set(legacy_RFL_INCLUDE_DIR "${legacy_RFL_SOURCE_DIR}/include")
add_subdirectory(${legacy_RFL_SOURCE_DIR}/source)

# Add core RFL library
set(core_RFL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/core")
set(core_RFL_INCLUDE_DIR "${core_RFL_SOURCE_DIR}")
add_subdirectory(${core_RFL_SOURCE_DIR})

# Add aliases for consumers
add_library(RFL::core ALIAS rfl_core)
add_library(RFL::rfl ALIAS rfl_core)
add_library(RFL::legacy ALIAS rfl_legacy)

# Add tests if building tests.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND NOT SKBUILD)
        add_subdirectory(${legacy_RFL_SOURCE_DIR}/tests)
endif()

if(SKBUILD)
    set(BUILD_PYTHON_BINDINGS_DEFAULT ON)
else()
    set(BUILD_PYTHON_BINDINGS_DEFAULT OFF)
endif()
option(BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" ${BUILD_PYTHON_BINDINGS_DEFAULT})
if(BUILD_PYTHON_BINDINGS)
    include(FetchContent)

    # Fetch carma (NumPy/Armadillo bridge)
    FetchContent_Declare(
        carma
        GIT_REPOSITORY https://github.com/RUrlus/carma.git
        GIT_TAG        v0.8.0
    )

    # Fetch pybind11
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.1
    )
    FetchContent_MakeAvailable(carma pybind11)

    add_subdirectory(python_bindings)
endif()