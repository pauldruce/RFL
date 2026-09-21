# Downstream Consumer Canary

This example demonstrates how an external C++ project consumes `RFL` via CMake `FetchContent`.

## CMake Integration

Add the following snippet to your downstream `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(MySimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
    rfl
    GIT_REPOSITORY https://github.com/pauldruce/RFL.git
    GIT_TAG        v0.3.0
)
FetchContent_MakeAvailable(rfl)

add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE RFL::core)
```

## Target Guarantees
* **Target Alias:** Linking `RFL::core` provides transitive include paths for all public headers.
* **Transitive Dependencies:** External consumers do not need manual include configurations for Armadillo (GSL is not required for `rfl_core`).
