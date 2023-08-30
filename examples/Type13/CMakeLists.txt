find_package(HDF5 REQUIRED)

add_executable(Type13Metropolis Type13Metropolis.cpp
        EigenvalueRecorder.hpp
        EigenvalueRecorder.cpp)

target_link_libraries(Type13Metropolis new_RFL HDF5::HDF5)
target_include_directories(Type13Metropolis PUBLIC ${new_RFL_INCUDE_DIR})
