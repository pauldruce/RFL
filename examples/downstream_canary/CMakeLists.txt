cmake_minimum_required(VERSION 3.20)
project(downstream_canary LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(TARGET RFL::core)
    # When built within RFL examples hierarchy
    add_executable(downstream_canary main.cpp)
    target_link_libraries(downstream_canary PRIVATE RFL::core)
else()
    # Standalone downstream consumer using CMake FetchContent
    include(FetchContent)
    if(DEFINED RFL_SOURCE_DIR)
        FetchContent_Declare(
            rfl
            SOURCE_DIR ${RFL_SOURCE_DIR}
        )
    else()
        FetchContent_Declare(
            rfl
            GIT_REPOSITORY https://github.com/pauldruce/RFL.git
            GIT_TAG main
        )
    endif()
    FetchContent_MakeAvailable(rfl)

    add_executable(downstream_canary main.cpp)
    target_link_libraries(downstream_canary PRIVATE RFL::core)
endif()
