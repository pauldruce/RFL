cc_binary(
    name = "Type13Metropolis",
    srcs = [
        "EigenvalueRecorder.cpp",
        "Type13Metropolis.cpp",
    ],
    hdrs = ["EigenvalueRecorder.hpp"],
    linker_flags = [
        "-lhdf5",
    ],
    deps = [
        "//src/RFL/new_RFL:rfl",
    ],
)