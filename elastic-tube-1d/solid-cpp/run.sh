#!/bin/bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build
  cmake --build build
fi

./build/SolidSolver ../precice-config.xml

close_log
