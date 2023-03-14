add_executable(MauroVsNew MauroVsNew.cpp)

target_include_directories(MauroVsNew PUBLIC ${new_RFL_INCUDE_DIR} ${RFL_INCLUDE_DIR})
target_link_libraries(MauroVsNew new_RFL RFL)