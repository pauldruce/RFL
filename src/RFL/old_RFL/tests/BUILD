package(cc = {
    "test_main": "///third_party/cc/gtest//:gmock",
})

cc_test(
    name = "unit_tests",
    srcs = [
        "action_tests.cpp",
        "clifford.test.cpp",
        "compilation_tests.cpp",
        "delta_tests.cpp",
        "derivative_tests.cpp",
        "geometry.test.cpp",
    ],
    deps = [
        "//src/RFL/old_RFL:rfl",
    ],
)