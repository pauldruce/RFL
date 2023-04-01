add_subdirectory(Type13)

add_executable(mauro_mmc_example
        mauro_thesis_example.cpp)

target_include_directories(mauro_mmc_example PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(mauro_mmc_example RFL)


add_executable(mauro_hmc_example
        mauro_hmc_example.cpp)

target_include_directories(mauro_hmc_example PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(mauro_hmc_example RFL)