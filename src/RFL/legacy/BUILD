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

# Export source files for the wheel build
filegroup(
    name = "wheel_srcs",
    srcs = glob(["**/*.cpp", "**/*.hpp", "**/*.h", "**/CMakeLists.txt"]),
    visibility = ["//src/RFL:rfl_wheel"],
)
