gentest(
    name = "unit_tests",
    data = [
        "//src/RFL:rfl_wheel",
    ] + glob(["unit/*.py"]),
    test_cmd = [
        "python3 -m venv _venv",
        "source _venv/bin/activate",
        "pip install --quiet --cache-dir /tmp/pip_cache numpy pytest",
        "pip install --quiet --cache-dir /tmp/pip_cache --no-index --find-links $(location //src/RFL:rfl_wheel) rfl",
        "pytest src/RFL/python_bindings/tests/unit",
    ],
    no_test_output = True,
    labels = ["unit"],
)

gentest(
    name = "integration_tests",
    data = [
        "//src/RFL:rfl_wheel",
    ] + glob(["integration/*.py"]),
    test_cmd = [
        "python3 -m venv _venv",
        "source _venv/bin/activate",
        "pip install --quiet --cache-dir /tmp/pip_cache numpy pytest",
        "pip install --quiet --cache-dir /tmp/pip_cache --no-index --find-links $(location //src/RFL:rfl_wheel) rfl",
        "pytest src/RFL/python_bindings/tests/integration",
    ],
    no_test_output = True,
    labels = ["integration"],
)
