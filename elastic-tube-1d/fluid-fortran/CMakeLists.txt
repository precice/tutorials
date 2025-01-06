cmake_minimum_required(VERSION 3.10)
project(ElasticTube LANGUAGES Fortran C)

# Set default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: ${CMAKE_BUILD_TYPE}")

# Find required packages
find_package(precice 3 REQUIRED CONFIG)
find_package(LAPACK REQUIRED)

# Add executable
add_executable(FluidSolver
  src/FluidComputeSolution.f90
  src/utilities.f90
  src/FluidSolver.f90
)

# Link libraries
target_link_libraries(FluidSolver PRIVATE precice::precice)
target_link_libraries(FluidSolver PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})

# ============================
# Add Compiler Flags for Debugging
# ============================

# Define a function to add compiler flags based on compiler ID
function(add_fortran_debug_flags target)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # GNU Fortran (gfortran) specific flags
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:Debug>:-fbounds-check>
      $<$<CONFIG:Debug>:-fbacktrace>
      $<$<CONFIG:Debug>:-fcheck=all>
      $<$<CONFIG:Debug>:-g>
      $<$<CONFIG:Debug>:-Wall>
      $<$<CONFIG:Debug>:-Wextra>
    )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # Intel Fortran (ifort) specific flags
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:Debug>:-check bounds>
      $<$<CONFIG:Debug>:-g>
      $<$<CONFIG:Debug>:-warn all>
    )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR
         CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
    # PGI, NAG, and other Fortran compilers
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:Debug>:-C>
      $<$<CONFIG:Debug>:-g>
      $<$<CONFIG:Debug>:-Wall>
    )
  else()
    message(WARNING "Compiler ${CMAKE_Fortran_COMPILER_ID} not explicitly supported for debug flags.")
    # Generic debug flags
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:Debug>:-g>
      $<$<CONFIG:Debug>:-Wall>
    )
  endif()
endfunction()

# Apply the debug flags to FluidSolver
add_fortran_debug_flags(FluidSolver)

# ============================
# Optionally, Add Linker Flags (if needed)
# ============================

# For example, to link LAPACK, you might need to specify additional flags
# However, typically, find_package(LAPACK) handles this
# If you have specific linker flags, you can add them as follows:
# target_link_libraries(FluidSolver PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})

# ============================
# Enable Fortran Module Generation (Optional)
# ============================

# If you have modules that need to be visible to other targets, set the appropriate properties
# For example:
# set_target_properties(FluidSolver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ============================
# Set CMake Policies (Optional)
# ============================

# To ensure compatibility with newer CMake versions
# cmake_policy(SET CMP0074 NEW)

# ============================
# End of CMakeLists.txt
# ============================
