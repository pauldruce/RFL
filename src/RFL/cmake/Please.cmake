# Set the output directory for vendors
set(VENDOR_OUTPUT_DIR "${GIT_PROJECT_ROOT}/plz-out/gen/vendors/")

# Add the vendor output directory to the CMake prefix path
list(INSERT CMAKE_PREFIX_PATH 0 "${VENDOR_OUTPUT_DIR}")

# Print the CMake prefix path for debugging
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")