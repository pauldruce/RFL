pybind11_add_module(rfl bindings.cpp)

# Find Python to execute script for Numpy headers
find_package(Python3 COMPONENTS Interpreter REQUIRED)
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Found NumPy Include Dir: ${NUMPY_INCLUDE_DIR}")
target_include_directories(rfl PRIVATE ${NUMPY_INCLUDE_DIR})

# Link to core rfl library, carma, and Armadillo
target_link_libraries(rfl PRIVATE rfl_core carma::carma)

# Pybind11 automatically requires C++11, but RFL needs C++17
target_compile_features(rfl PRIVATE cxx_std_17)

if(DEFINED SKBUILD_PROJECT_VERSION_FULL)
    target_compile_definitions(rfl PRIVATE RFL_VERSION_STRING="${SKBUILD_PROJECT_VERSION_FULL}")
elseif(DEFINED PROJECT_VERSION)
    target_compile_definitions(rfl PRIVATE RFL_VERSION_STRING="${PROJECT_VERSION}")
endif()

# Set the install location so pip can bundle it correctly
install(TARGETS rfl DESTINATION .)
