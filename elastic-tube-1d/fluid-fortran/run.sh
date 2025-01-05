#!/usr/bin/env bash
set -e -u

if [ ! -d build ]; then
  mkdir build
  cd build
  cmake ..
  cmake --build .
  cd ..
fi

./build/FluidSolver ../precice-config.xml
