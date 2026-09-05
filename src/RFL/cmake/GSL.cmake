# The purpose of this CMake module is to find the GNU Scientific Library (GSL).
# It first checks if GSL has already been defined as a target.
# If not, it uses CMake's standard find_package(GSL REQUIRED).

if(NOT TARGET GSL::gsl)
    find_package(GSL REQUIRED)
endif()

if(TARGET GSL::gsl)
    set(GSL_LIBRARIES GSL::gsl)
    # cspell:disable-next-line
    if(TARGET GSL::gslcblas)
        # cspell:disable-next-line
        list(APPEND GSL_LIBRARIES GSL::gslcblas)
    endif()
    set(GSL_FOUND TRUE)
elseif(GSL_FOUND)
    # Ensure modern imported target exists even on legacy find modules
    if(NOT TARGET GSL::gsl)
        add_library(GSL::gsl INTERFACE IMPORTED)
        if(GSL_INCLUDE_DIRS)
            target_include_directories(GSL::gsl INTERFACE ${GSL_INCLUDE_DIRS})
        endif()
        if(GSL_LIBRARIES)
            target_link_libraries(GSL::gsl INTERFACE ${GSL_LIBRARIES})
        endif()
    endif()
    set(GSL_LIBRARIES GSL::gsl)
else()
    message(FATAL_ERROR "GNU Scientific Library (GSL) not found. Please install libgsl-dev (Linux) or brew install gsl (macOS).")
endif()

if(DEFINED GSL_INCLUDE_DIRS)
    include_directories(${GSL_INCLUDE_DIRS})
endif()