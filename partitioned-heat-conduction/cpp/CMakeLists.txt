cmake_minimum_required (VERSION 3.10)

project(Heat LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)

find_package(precice REQUIRED CONFIG)

add_executable(heat heat.cpp)
target_link_libraries(heat PRIVATE precice::precice)
