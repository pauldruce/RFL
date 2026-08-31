add_library(rfl_core
        IDiracOperator.hpp
        IDiracOperatorDerivatives.cpp
        IDiracOperatorDerivatives.hpp)

set(rfl_core_HEADERS
        BarrettGlaser/Action.hpp
        BarrettGlaser/Hamiltonian.hpp
        BarrettGlaser/Metropolis.hpp
        Clifford.hpp
        DiracOperator.hpp
        GslRng.hpp
        IAction.hpp
        IAlgorithm.hpp
        IRng.hpp
        Simulation.hpp)

target_sources(
        rfl_core
        PRIVATE
        BarrettGlaser/Action.cpp
        BarrettGlaser/Hamiltonian.cpp
        BarrettGlaser/Metropolis.cpp
        Clifford.cpp
        DiracOperator.cpp
        Simulation.cpp
        PUBLIC ${rfl_core_HEADERS})

target_include_directories(rfl_core PUBLIC ${core_RFL_INCLUDE_DIR})
target_link_libraries(rfl_core PUBLIC ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})

# Install static library files
install(
        TARGETS rfl_core
        DESTINATION lib/rfl/bin/${CMAKE_BUILD_TYPE}
)

# Install include/header files.
install(
        DIRECTORY ${core_RFL_INCLUDE_DIR}/
        DESTINATION lib/rfl/include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "tests" EXCLUDE
)

if(NOT SKBUILD)
    add_subdirectory(tests)
endif()
