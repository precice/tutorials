#!/bin/sh
set -e -u

. ../../tools/log.sh

if [ ! -d build ]; then
  log mkdir build
  log cmake -S . -B build
  log cmake --build build
fi

log ./build/FluidSolver ../precice-config.xml

close_log
