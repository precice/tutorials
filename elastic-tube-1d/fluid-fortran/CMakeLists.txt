cmake_minimum_required(VERSION 3.10)
project(ElasticTube_Fluid_Fortran LANGUAGES C Fortran DESCRIPTION "preCICE Fortran Fluid Solver")

# If no build type is specified, default to Debug:
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: ${CMAKE_BUILD_TYPE}")

# Find preCICE library with Fortran bindings
find_package(precice 3 REQUIRED CONFIG)

# Find LAPACK for linear solves
# Uncomment if LAPACK is needed
# find_package(LAPACK REQUIRED)

# Add the executable
add_executable(FluidSolver
# Uncomment these if additional source files are required
#   src/FluidComputeSolution.f90
#   src/utilities.f90
  src/FluidSolver.f90
)

# Link preCICE library
target_link_libraries(FluidSolver PRIVATE precice::precice)

# Additional linking for LAPACK if needed
# target_link_libraries(FluidSolver PRIVATE LAPACK::LAPACK)
