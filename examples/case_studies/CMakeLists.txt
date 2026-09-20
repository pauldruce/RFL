find_package(HDF5 QUIET)
if(HDF5_FOUND)
    add_subdirectory(type_13_simulation)
else()
    message(STATUS "HDF5 not found: skipping type_13_simulation case study")
endif()

add_subdirectory(mauro_thesis_mmc)
add_subdirectory(hmc_tuning)
