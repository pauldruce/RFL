cmake_minimum_required(VERSION 3.7...3.23)

if(${CMAKE_VERSION} VERSION_LESS 3.23 )
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
        cmake_policy(VERSION 3.23)
endif()

project(RFL
        VERSION 1.0
        DESCRIPTION
        "RFL - A library for Markov Chain Monte Carlo simulations of Finite Random Non-commutative geometries."
        LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STADARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#     include(CTest)
endif()

# Get DOCTEST
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
        doctest
        PREFIX ${CMAKE_BINARY_DIR}/doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        TIMEOUT 10
        UPDATE_COMMAND ${GIT_EXECUTABLE} pull
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON)

# Expose required variable (DOCTEST_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(doctest source_dir)
set(DOCTEST_INCLUDE_DIR
        ${source_dir}/doctest
        CACHE INTERNAL "Path to include folder for doctest")
include_directories(DOCTEST_INCLUDE_DIR)

# Build library
set(INCLUDE_DIR ${RFL_SOURCE_DIR}/include)
add_subdirectory(source)
# Add tests if building tests.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(tests)
endif()
