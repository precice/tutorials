cmake_minimum_required(VERSION 3.10)
project(ElasticTube_Solid_Fortran LANGUAGES C Fortran DESCRIPTION "preCICE Fortran Solid Solver")


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: ${CMAKE_BUILD_TYPE}")


find_package(precice REQUIRED CONFIG)


add_executable(SolidSolver
  src/SolidSolver.f90
)


target_link_libraries(SolidSolver PRIVATE precice::precice)
