python_test(
    name = "unit_tests",
    srcs = glob(["unit/test_*.py"]),
    deps = [
        "//src/RFL:rfl_wheel",
        "//third_party/python:pytest",
        "//third_party/python:numpy"
    ],
    labels = ["unit"],
)

python_test(
    name = "integration_tests",
    srcs = glob(["integration/test_*.py"]),
    deps = [
        "//src/RFL:rfl_wheel",
        "//third_party/python:pytest",
        "//third_party/python:numpy"
    ],
    labels = ["integration"],
)
