#copy all source code to build directory to make available to pybind
add_custom_target(copy_micro_heat_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# add a symlink for each input file
add_input_file_links()
    
#build using pybind
add_subdirectory(pybind11)
pybind11_add_module(micro_sim micro_sim.cpp)
target_link_libraries(micro_sim PRIVATE dumux-precice)