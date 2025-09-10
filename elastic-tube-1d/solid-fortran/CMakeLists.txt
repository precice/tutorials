cmake_minimum_required(VERSION 3.10)

project(ElasticTube LANGUAGES Fortran C)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -fbounds-check")
endif()

find_package(precice 3 REQUIRED CONFIG)

add_executable(SolidSolver
  src/SolidComputeSolution.f90
  src/SolidSolver.f90
)

target_link_libraries(SolidSolver PRIVATE precice::precice)
