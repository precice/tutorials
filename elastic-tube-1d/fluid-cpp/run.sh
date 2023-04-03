#!/bin/sh
set -e -u

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build
  cmake --build build
fi

./build/FluidSolver ../precice-config.xml
