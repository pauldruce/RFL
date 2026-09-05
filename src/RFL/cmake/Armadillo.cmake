# The purpose of this CMake module is to find or fetch the Armadillo library.
# It first checks if Armadillo has already been defined as a target.
# If not, and unless RFL_FETCH_ARMADILLO is enabled, it attempts find_package(Armadillo QUIET).
# If find_package fails or RFL_FETCH_ARMADILLO is ON, it falls back to FetchContent
# from the upstream GitLab repository.

if(NOT TARGET Armadillo::Armadillo AND NOT TARGET armadillo)
    if(NOT RFL_FETCH_ARMADILLO)
        message(STATUS "Searching for system Armadillo library...")
        find_package(Armadillo QUIET)
    endif()

    if(NOT Armadillo_FOUND AND NOT ARMADILLO_FOUND)
        message(STATUS "Armadillo not found or RFL_FETCH_ARMADILLO is enabled; fetching via FetchContent...")
        include(FetchContent)

        if(NOT DEFINED RFL_ARMADILLO_GIT_TAG)
            if(DEFINED ARMADILLO_GIT_TAG)
                set(RFL_ARMADILLO_GIT_TAG "${ARMADILLO_GIT_TAG}")
            else()
                set(RFL_ARMADILLO_GIT_TAG "14.2.x")
            endif()
        endif()

        set(BUILD_SMOKE_TEST OFF CACHE BOOL "Disable Armadillo smoke test" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static Armadillo library" FORCE)

        if(APPLE)
            set(ALLOW_OPENBLAS_MACOS ON CACHE BOOL "Allow detection of OpenBLAS on macOS" FORCE)
        endif()

        FetchContent_Declare(
            armadillo
            GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
            GIT_TAG        ${RFL_ARMADILLO_GIT_TAG}
        )
        FetchContent_MakeAvailable(armadillo)
    endif()
endif()

if(TARGET armadillo)
    if(NOT TARGET Armadillo::Armadillo)
        add_library(Armadillo::Armadillo ALIAS armadillo)
    endif()
    set(ARMADILLO_LIBRARIES Armadillo::Armadillo)
    set(ARMADILLO_FOUND TRUE)
elseif(TARGET Armadillo::Armadillo)
    set(ARMADILLO_LIBRARIES Armadillo::Armadillo)
    set(ARMADILLO_FOUND TRUE)
elseif(ARMADILLO_FOUND OR Armadillo_FOUND)
    set(ARMADILLO_FOUND TRUE)
    if(NOT TARGET Armadillo::Armadillo)
        add_library(Armadillo::Armadillo INTERFACE IMPORTED)
        if(ARMADILLO_INCLUDE_DIRS)
            target_include_directories(Armadillo::Armadillo INTERFACE ${ARMADILLO_INCLUDE_DIRS})
        endif()
        if(ARMADILLO_LIBRARIES)
            target_link_libraries(Armadillo::Armadillo INTERFACE ${ARMADILLO_LIBRARIES})
        endif()
    endif()
    set(ARMADILLO_LIBRARIES Armadillo::Armadillo)
else()
    message(FATAL_ERROR "Armadillo library could not be found or fetched.")
endif()

if(DEFINED ARMADILLO_INCLUDE_DIRS)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
endif()
