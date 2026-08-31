find_package(HDF5 REQUIRED)

add_executable(Type13Metropolis Type13Metropolis.cpp
        EigenvalueRecorder.hpp
        EigenvalueRecorder.cpp)

target_link_libraries(Type13Metropolis PRIVATE rfl_core HDF5::HDF5)
target_include_directories(Type13Metropolis PUBLIC ${core_RFL_INCLUDE_DIR})
