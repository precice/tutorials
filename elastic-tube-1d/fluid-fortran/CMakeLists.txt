cmake_minimum_required(VERSION 3.10)
project(ElasticTube LANGUAGES Fortran C)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: ${CMAKE_BUILD_TYPE}")

find_package(precice 3 REQUIRED CONFIG)
find_package(LAPACK REQUIRED)

add_executable(FluidSolver
  src/FluidComputeSolution.f90
  src/utilities.f90
  src/FluidSolver.f90
)

target_link_libraries(FluidSolver PRIVATE precice::precice)
target_link_libraries(FluidSolver PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})

