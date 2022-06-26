cmake_minimum_required(VERSION 3.7...3.23)

if (${CMAKE_VERSION} VERSION_LESS 3.23)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else ()
    cmake_policy(VERSION 3.23)
endif ()

project(
        RFL
        VERSION 1.0
        DESCRIPTION
        "RFL - A library for Markov Chain Monte Carlo simulations of Finite Random Non-commutative geometries."
        LANGUAGES CXX)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif ()

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# GSL is included with Armadillo above, but we need more?
find_package(GSL REQUIRED)

### Get DOCTEST
include(FetchContent)
FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG v2.4.8)
FetchContent_MakeAvailable(doctest)
set(DOCTEST_INCLUDE_DIR
        "${doctest_SOURCE_DIR}/"
        CACHE INTERNAL "Path to include folder for doctest")

include_directories(DOCTEST_INCLUDE_DIR)

### Get GoogleTest
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)



### Build library
set(INCLUDE_DIR ${RFL_SOURCE_DIR}/include)
add_subdirectory(source)

### Add tests if building tests.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif ()
