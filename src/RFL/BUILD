genrule(
    name = "rfl_wheel",
    srcs = [
        "pyproject.toml",
        "CMakeLists.txt",
        "python_bindings/bindings.cpp",
        "python_bindings/CMakeLists.txt",
        # Sub-package source files (exported via wheel_srcs filegroups)
        "//src/RFL/new_RFL:wheel_srcs",
        "//src/RFL/old_RFL:wheel_srcs",
    ] + glob(["cmake/*.cmake"]),
    outs = ["rfl.whl"],
    cmd = [
        # Set up build tools in a venv
        "python3 -m venv _venv",
        "source _venv/bin/activate",
        "pip install --quiet --cache-dir /tmp/pip_cache --upgrade pip",
        "pip install --quiet --cache-dir /tmp/pip_cache build scikit-build-core pybind11 cmake ninja",
        # Build the wheel (srcs are linked into TMP_DIR preserving paths, command runs in TMP_DIR)
        "python -m build --wheel --outdir _wheeldir $PKG_DIR",
        # Copy the built wheel to the output location
        "cp _wheeldir/*.whl $OUTS",
    ],
    timeout = 1800, # 30 minutes to allow for C++ compilation
    visibility = ["PUBLIC"],
)
