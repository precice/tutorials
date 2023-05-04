add_executable(test_macro_heat main.cc)
target_compile_definitions(test_macro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_macro_heat PRIVATE dumux-precice)
target_compile_definitions(test_macro_heat PUBLIC TYPETAG=OnePNIConductionCCTpfa)

# add a symlink for each input file
add_input_file_links()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/clean-example.sh
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
