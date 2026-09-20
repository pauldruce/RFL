add_executable(hmc_tuning hmc_tuning.cpp)
target_include_directories(hmc_tuning PUBLIC ${legacy_RFL_INCLUDE_DIR})
target_link_libraries(hmc_tuning PRIVATE rfl_legacy)
