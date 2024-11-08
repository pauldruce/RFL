cc_library(
    name = "gtest",
    srcs = glob(
        ["googletest/src/*.cc"],
        exclude = ["googletest/src/gtest-all.cc"],
    ),
    hdrs = glob(["googletest/include/**.h"]),
    includes = [
        "googletest",
        "googletest/include",
    ],
    linker_flags = ["-lpthread"],
    private_hdrs = glob(["googletest/src/*.h"]),
    visibility = ["PUBLIC"],
)

cc_library(
    name = "gmock",
    srcs = glob(
        include = [
            "googlemock/src/*.cc",
        ],
        exclude = [
            "googlemock/src/gmock-all.cc",
            "googlemock/src/gmock_main.cc",
        ],
    ),
    hdrs = glob([
        "googlemock/include/gmock/**.h",
    ]),
    private_hdrs = glob([
        "googlemock/include/gmock/**/*.h",
    ]),
    includes = [
        "googlemock",
        "googlemock/include",
    ],
    deps = [":gtest"],
    linker_flags = ["-lpthread"],
    visibility = ["PUBLIC"],
)
