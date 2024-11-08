package(cc = {
    "test_main": "///third_party/cc/gtest//:gmock",
})

cc_test(
    name = "unit_tests",
    srcs = [
        "./tAction.cpp",
        "./tClifford.cpp",
        "./tDiracOperator.cpp",
        "./tGslRng.cpp",
        "./tHamiltonian.cpp",
        "./tMetropolis.cpp",
        "./tSimulation.cpp",
    ],
    deps = [
        "//src/RFL/new_RFL:rfl",
    ],
)
