cmake_minimum_required(VERSION 3.7...3.23)

if(${CMAKE_VERSION} VERSION_LESS 3.23)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.23)
endif()

project(
        RFL
        VERSION 0.1
        DESCRIPTION
        "Random Fuzzy Library (RFL) - A library for Markov Chain Monte Carlo simulations of Finite Random Non-commutative geometries."
        LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  enable_testing()

  add_compile_options(-Wall -Wextra -Wpedantic )
  add_compile_options(-fsanitize=address )
  add_link_options(-fsanitize=address )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(Doxygen)
find_package(Doxygen
        REQUIRED dot)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_HTML_OUTPUT
        ${PROJECT_SOURCE_DIR}/RFL_docs)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE
        "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
# Settings for Inheritance Diagram etc.
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_CLASS_DIAGRAMS YES)
set(DOXYGEN_HIDE_UNDOC_RELATIONS NO)
set(DOXYGEN_CALL_GRAPH YES)
set(DOXYGEN_CALLER_GRAPH YES)
set(DOXYGEN_HAVE_DOT             YES)
set(DOXYGEN_CLASS_GRAPH          YES)
set(DOXYGEN_COLLABORATION_GRAPH  NO)
set(DOXYGEN_UML_LOOK             YES)
set(DOXYGEN_UML_LIMIT_NUM_FIELDS 50)
set(DOXYGEN_TEMPLATE_RELATIONS   YES)
set(DOXYGEN_DOT_GRAPH_MAX_NODES  100)
set(DOXYGEN_MAX_DOT_GRAPH_DEPTH  0)
set(DOXYGEN_DOT_TRANSPARENT      YES)
set(DOXYGEN_INTERACTIVE_SVG  YES)
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_EXCLUDE_PATTERNS */tests/*)

doxygen_add_docs(doxygen
        ${PROJECT_SOURCE_DIR}/README.md ${PROJECT_SOURCE_DIR}/RFL_source/new_source ${PROJECT_SOURCE_DIR}/examples
        COMMENT "Generate documentation for entire project"
        )

# GSL is included with Armadillo above, but we need more.
find_package(Armadillo REQUIRED)
find_package(GSL REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
add_subdirectory(vendors)

add_subdirectory(RFL_source)
add_subdirectory(examples)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


# AN EXPERIMENT
add_subdirectory(playground)
