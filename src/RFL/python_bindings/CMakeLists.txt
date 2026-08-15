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

# Link to new_RFL, carma, and Armadillo (which new_RFL uses)
target_link_libraries(rfl PRIVATE new_RFL carma::carma)

# Pybind11 automatically requires C++11, but RFL needs C++17
target_compile_features(rfl PRIVATE cxx_std_17)

# Set the install location so pip can bundle it correctly
install(TARGETS rfl DESTINATION .)
