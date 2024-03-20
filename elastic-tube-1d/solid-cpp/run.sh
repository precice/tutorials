#!/bin/sh
set -e -u

. ../../tools/log.sh

if [ ! -d build ]; then
  mkdir build
  log cmake -S . -B build
  log cmake --build build
fi

log ./build/SolidSolver ../precice-config.xml

close_log
