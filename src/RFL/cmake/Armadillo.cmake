# The purpose of this Cmake module is to find the Armadillo library.
# This library may be installed by the system/user in the standard install
# location, or it might be installed to a custom location by the build system
# Please.
set(LOCAL_ARMADILLO_INSTALL_DIR "${GIT_PROJECT_ROOT}/plz-out/gen/vendors/armadillo-12.6.6/install")

# Create local variable which stores the files to search for
message(DEBUG "Checking for local Armadillo library at ${LOCAL_ARMADILLO_INSTALL_DIR}")

if(EXISTS "${LOCAL_ARMADILLO_INSTALL_DIR}/lib/libarmadillo${CMAKE_SHARED_LIBRARY_SUFFIX}"
    OR EXISTS "${LOCAL_ARMADILLO_INSTALL_DIR}/lib/libarmadillo${CMAKE_STATIC_LIBRARY_SUFFIX}")
    message(STATUS "Using local Armadillo library at ${LOCAL_ARMADILLO_INSTALL_DIR}")
    set(ARMADILLO_LIBRARIES "${LOCAL_ARMADILLO_INSTALL_DIR}/lib/libarmadillo${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(ARMADILLO_INCLUDE_DIRS "${LOCAL_ARMADILLO_INSTALL_DIR}/include")

    set(ARMADILLO_FOUND TRUE)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    link_directories(${ARMADILLO_LIBRARIES})
else()
    message(STATUS "Using system Armadillo library")
    find_package(Armadillo REQUIRED)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
endif()

# If the GSL package is not found, stop with an error
if(NOT ARMADILLO_FOUND)
    message(FATAL_ERROR "Armadillo library not found")
endif()
