#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh || true
exec > >(tee --append "$LOGFILE") 2>&1 || true

if [ ! -f src/precice.f90 ]; then
  echo "Fetching precice.f90 (Module for Fortran bindings of preCICE)..."
  curl -o src/precice.f90 https://raw.githubusercontent.com/precice/fortran-module/master/precice.f90
fi

if [ ! -d build ]; then
  mkdir build
  cd build
  cmake ..
  cmake --build .
  cd ..
fi

./build/SolidSolver ../precice-config.xml

close_log
