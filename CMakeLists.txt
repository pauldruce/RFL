cmake_minimum_required(VERSION 3.7...3.23)

if (${CMAKE_VERSION} VERSION_LESS 3.23)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else ()
    cmake_policy(VERSION 3.23)
endif ()

project(
        RFL
        VERSION 0.1
        DESCRIPTION
        "RFL = Random Fuzzy Library - A library for Markov Chain Monte Carlo simulations of Finite Random Non-commutative geometries."
        LANGUAGES CXX)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif ()


# GSL is included with Armadillo above, but we need more.
find_package(Armadillo REQUIRED)
find_package(GSL REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
add_subdirectory(vendors)

add_subdirectory(RFL_source)
add_subdirectory(examples)
