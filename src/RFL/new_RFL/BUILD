cc_library(
    name = "rfl",
    srcs = [
        "BarrettGlaser/Action.cpp",
        "BarrettGlaser/Hamiltonian.cpp",
        "BarrettGlaser/Metropolis.cpp",
        "Clifford.cpp",
        "DiracOperator.cpp",
        "IDiracOperatorDerivatives.cpp",
        "Simulation.cpp",
    ],
    hdrs = [
        "BarrettGlaser/Action.hpp",
        "BarrettGlaser/Hamiltonian.hpp",
        "BarrettGlaser/Metropolis.hpp",
        "Clifford.hpp",
        "DiracOperator.hpp",
        "GslRng.hpp",
        "IAction.hpp",
        "IAlgorithm.hpp",
        "IDiracOperator.hpp",
        "IDiracOperatorDerivatives.hpp",
        "IRng.hpp",
        "Simulation.hpp",
    ],
    includes = [
        "./",  # TODO: This is not ideal, I should keep my header files in an include directory which matches what the output will be.
    ],
    linker_flags = [
        "-larmadillo",  # Depends on armadillo
        "-lgsl",  # Depends on GSL
    ],
    visibility = ["PUBLIC"],
)