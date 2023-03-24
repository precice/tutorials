#!/bin/sh
set -e -u

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build
fi

cmake --build build

./build/SolidSolver ../precice-config.xml
