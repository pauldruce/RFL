add_library(new_RFL
        IDiracOperator.hpp
        IDiracOperatorDerivatives.cpp
        IDiracOperatorDerivatives.hpp)

set(new_RFL_HEADERS
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
        new_RFL
        PRIVATE
        BarrettGlaser/Action.cpp
        BarrettGlaser/Hamiltonian.cpp
        BarrettGlaser/Metropolis.cpp
        Clifford.cpp
        DiracOperator.cpp
        Simulation.cpp
        PUBLIC ${new_RFL_HEADERS})

target_include_directories(new_RFL PUBLIC ${new_RFL_INCLUDE_DIR})
target_link_libraries(new_RFL ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})

# Install static library files
install(
        TARGETS new_RFL
        DESTINATION ${CMAKE_SOURCE_DIR}/lib/new_RFL/bin/${CMAKE_BUILD_TYPE}
)

# Install include/header files.
install(
        DIRECTORY ${new_RFL_INCLUDE_DIR}/
        DESTINATION ${CMAKE_SOURCE_DIR}/lib/new_RFL/include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "tests" EXCLUDE
)

add_subdirectory(tests)
