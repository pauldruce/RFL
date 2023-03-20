add_library(new_RFL )

set(new_RFL_HEADERS
        Action.hpp
        Clifford.hpp
        DiracOperator.hpp
        GslRng.hpp
        Hamiltonian.hpp
        IAction.hpp
        IAlgorithm.hpp
        IRng.hpp
        Metropolis.hpp
        Simulation.hpp)

target_sources(
        new_RFL
        PRIVATE
        Action.cpp
        Clifford.cpp
        DiracOperator.cpp
        Hamiltonian.cpp
        Metropolis.cpp
        Simulation.cpp
        PUBLIC ${new_RFL_HEADERS})

set_target_properties(new_RFL PROPERTIES PUBLIC_HEADER "${new_RFL_HEADERS}")
target_include_directories(new_RFL PUBLIC ${new_RFL_INCLUDE_DIR})
target_link_libraries(new_RFL ${ARMADILLO_LIBRARIES} ${GSL_LIBRARIES})

install(
        TARGETS new_RFL
        DESTINATION ${CMAKE_SOURCE_DIR}/lib/${CMAKE_CONFIG}/new_RFL/bin
        PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/lib/new_RFL/include
)

add_subdirectory(tests)
