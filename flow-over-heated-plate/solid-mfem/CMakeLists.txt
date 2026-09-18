cmake_minimum_required(VERSION 3.13)

project(solid-mfem LANGUAGES CXX)

find_package(MFEM 4.10 REQUIRED)
find_package(precice 3 REQUIRED)

add_executable(solid solid.cpp)
target_link_libraries(solid PRIVATE mfem precice::precice)
target_compile_features(solid PUBLIC cxx_std_17)
