# We require version CMake version 3.1 to prevent issues
# with dune_enable_all_packages and older CMake versions.
cmake_minimum_required(VERSION 3.1)
project(dumux-phasefield CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find_package(pybind11 CONFIG REQUIRED)
#pybind11_add_module(dumux-phasefield test/microsim/micro_sim.cpp)
#add_subdirectory(pybind11)

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages()

add_subdirectory(doc)
add_subdirectory(dumux)
add_subdirectory(test)
add_subdirectory(cmake/modules)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
