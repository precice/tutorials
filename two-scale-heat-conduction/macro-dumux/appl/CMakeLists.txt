add_executable(macro_dumux main.cc)
target_compile_definitions(macro_dumux PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(macro_dumux PRIVATE dumux-precice)
target_compile_definitions(macro_dumux PUBLIC TYPETAG=OnePNIConductionCCTpfa)

# add a symlink for each input file
add_input_file_links()
