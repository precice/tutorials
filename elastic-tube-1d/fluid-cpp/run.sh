#!/bin/bash
set -e -u

. ../../tools/log.sh

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build
  cmake --build build
fi

./build/FluidSolver ../precice-config.xml

close_log
