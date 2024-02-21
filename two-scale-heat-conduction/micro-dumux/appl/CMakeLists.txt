execute_process(COMMAND pip show pybind11
                COMMAND grep Location 
                COMMAND cut -d " " -f2
                OUTPUT_VARIABLE SITE_PACKAGES_DIR 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set(pybind11_DIR "${SITE_PACKAGES_DIR}/pybind11/share/cmake/pybind11/")

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(micro_sim micro_sim.cpp)

#copy all source code to build directory to make available to pybind
add_custom_target(copy_micro_heat_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
    
# add a symlink for each input file
add_input_file_links()

