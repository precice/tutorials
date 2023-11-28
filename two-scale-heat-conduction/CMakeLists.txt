project(coupled-heat-conduction-dumux-macro CXX)

cmake_minimum_required(VERSION 3.22)

add_executable(macro_dumux EXCLUDE_FROM_ALL main.cc)

find_package(dumux-precice REQUIRED)

target_link_libraries(macro_dumux PRIVATE dumux-precice)

message(STATUS "Making dune available")
function(clone_dune NAME)
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/_dune/dune-${NAME})
        execute_process(COMMAND git clone --depth=1 https://gitlab.dune-project.org/core/dune-${NAME}.git ${CMAKE_BINARY_DIR}/_dune/dune-${NAME})
    endif()
endfunction()

clone_dune(common)
clone_dune(geometry)
clone_dune(grid)
file(
    WRITE ${CMAKE_BINARY_DIR}/_dune/cmake.opts
    "CMAKE_FLAGS=\"-DDUNE_ENABLE_PYTHONBINDINGS=0  # do you need bindings?
                   -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/_dune_install\"”  # this depends on if you want to “install” or not...
)
execute_process(
    COMMAND ./dune-common/bin/dunecontrol --opts=cmake.opts all
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/_dune
)
# depending on if you want to install, you may need the following:
execute_process(
    COMMAND ./dune-common/bin/dunecontrol --opts=cmake.opts make install
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/_dune
)
