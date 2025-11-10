add_executable(free_flow_dumux main.cc)
target_compile_definitions(free_flow_dumux PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(free_flow_dumux PRIVATE dumux-precice)

add_input_file_links()