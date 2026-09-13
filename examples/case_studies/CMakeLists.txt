add_subdirectory(type_13_simulation)

add_executable(mauro_thesis_mmc mauro_thesis_mmc.cpp)
target_include_directories(mauro_thesis_mmc PUBLIC ${legacy_RFL_INCLUDE_DIR})
target_link_libraries(mauro_thesis_mmc PRIVATE rfl_legacy)

add_executable(hmc_tuning hmc_tuning.cpp)
target_include_directories(hmc_tuning PUBLIC ${legacy_RFL_INCLUDE_DIR})
target_link_libraries(hmc_tuning PRIVATE rfl_legacy)
