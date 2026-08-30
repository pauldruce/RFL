add_executable(main main.cpp)
target_include_directories(main PUBLIC ${RFL_INCLUDE_DIR})
target_link_libraries(main RFL)
