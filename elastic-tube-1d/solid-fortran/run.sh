#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh || true
exec > >(tee --append "$LOGFILE") 2>&1 || true

if [ ! -d build ]; then
  mkdir build
  cd build
  cmake ..
  cmake --build .
  cd ..
fi

./build/SolidSolver ../precice-config.xml

close_log
