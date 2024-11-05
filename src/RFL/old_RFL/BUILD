cc_library(
    name = "rfl",
    srcs = glob(["source/**/*.cpp"]),
    hdrs = glob(["include/**/*.hpp"]),
    includes = ["include"],
    linker_flags = [
        "-larmadillo",  # Depends on armadillo
        "-lgsl",  # Depends on GSL
    ],
    visibility = ["//src/RFL/..."],
)