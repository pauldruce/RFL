add_executable(main main.cpp)
target_include_directories(main PUBLIC ${core_RFL_INCLUDE_DIR})
target_link_libraries(main PRIVATE rfl_core)
