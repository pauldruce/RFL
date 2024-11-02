# Set the path to the local GSL installation
set(LOCAL_GSL_PATH "${GIT_PROJECT_ROOT}/plz-out/gen/vendors/gsl-2.8/install")

# Check if the GSL library is available in the local installation path
message(DEBUG "Checking for local GSL library at ${LOCAL_GSL_PATH}")

if(EXISTS "${LOCAL_GSL_PATH}/lib/libgsl${CMAKE_STATIC_LIBRARY_SUFFIX}"
    OR EXISTS "${LOCAL_GSL_PATH}/lib/libgsl${CMAKE_SHARED_LIBRARY_SUFFIX}")
    message(STATUS "Using local GSL library at ${LOCAL_GSL_PATH}")

    # Set GSL include and library paths
    set(GSL_INCLUDE_DIR "${LOCAL_GSL_PATH}/include")
    set(GSL_LIBRARIES "${LOCAL_GSL_PATH}/lib/libgsl${CMAKE_STATIC_LIBRARY_SUFFIX}")

    # Set GSL_FOUND to true as the library is locally found
    set(GSL_FOUND TRUE)

    # Add the GSL include directory to the include path
    include_directories(${GSL_INCLUDE_DIR})

    # Link the local GSL library
    link_directories(${GSL_LIBRARIES})
else()
    # Fallback to system GSL library using find_package
    message(STATUS "Using system GSL library")
    find_package(GSL REQUIRED)
    include_directories(${GSL_INCLUDE_DIRS})
endif()

# If the GSL package is not found, stop with an error
if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL library not found")
endif()