genrule(
    name = "build_armadillo",
    srcs = glob(["armadillo-12.6.6/**/*"]),
    cmd = [
        "cmake -B ./build -DALLOW_OPENBLAS_MACOS=ON vendors/armadillo-12.6.6",
        "cmake --build ./build --target all -j 4",
        "cmake --install ./build --prefix armadillo-12.6.6/install",
    ],
    outs = ["armadillo-12.6.6/"],
    visibility = ["//visibility:public"],
    tools = ["cmake"],
)


genrule(
    name = "build_gsl",
    srcs = glob(["gsl-2.8/**/*"]),
    tools = ["make"],
    cmd = [
        "cd vendors/gsl-2.8",
        "./configure --prefix=${HOME}/gsl-2.8/install",
        "make",
        "make install ",
    ],
    outs= ["gsl-2.8/"],
)