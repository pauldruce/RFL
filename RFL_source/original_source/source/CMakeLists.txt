add_library(RFL STATIC)

set(RFL_INCLUDES ${RFL_INCLUDE_DIR}/Cliff.hpp ${RFL_INCLUDE_DIR}/Geom24.hpp)
target_sources(
        RFL
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
        PUBLIC ${RFL_INCLUDES})

set_target_properties(RFL PROPERTIES PUBLIC_HEADER "${RFL_INCLUDES}")
target_include_directories(RFL PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(RFL ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})

# Install static library and header files
install(TARGETS RFL
        DESTINATION ${CMAKE_SOURCE_DIR}/lib/RFL/bin/${CMAKE_BUILD_TYPE}
        PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/lib/RFL/include
)
