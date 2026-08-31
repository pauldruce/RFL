add_library(rfl_legacy STATIC)

set(rfl_legacy_INCLUDES ${legacy_RFL_INCLUDE_DIR}/Cliff.hpp ${legacy_RFL_INCLUDE_DIR}/Geom24.hpp)
target_sources(
        rfl_legacy
        PRIVATE action.cpp
        cliff.cpp
        derivative.cpp
        hamiltonian.cpp
        HMC.cpp
        HMC_core.cpp
        leapfrog.cpp
        metropolis.cpp
        misc.cpp
        MMC.cpp
        MMC_core.cpp
        PUBLIC ${rfl_legacy_INCLUDES})

set_target_properties(rfl_legacy PROPERTIES PUBLIC_HEADER "${rfl_legacy_INCLUDES}")
target_include_directories(rfl_legacy PUBLIC ${legacy_RFL_INCLUDE_DIR})
target_link_libraries(rfl_legacy ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})

# Install static library and header files
install(TARGETS rfl_legacy
        DESTINATION lib/rfl_legacy/bin/${CMAKE_BUILD_TYPE}
        PUBLIC_HEADER DESTINATION lib/rfl_legacy/include
)
