add_executable(porous_media_dumux main.cc)
target_compile_definitions(porous_media_dumux PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(porous_media_dumux PRIVATE dumux-precice)

add_input_file_links()