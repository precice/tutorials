#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ ! -f src/precice.f90 ]; then
  # Get the preCICE Fortran module. Switch the branch with ./run.sh <branch>.
  DEFAULT_BRANCH="main"
  echo "Fetching precice.f90 from  https://github.com/precice/fortran-module/tree/${1:-$DEFAULT_BRANCH}..."
  curl -o src/precice.f90 "https://raw.githubusercontent.com/precice/fortran-module/${1:-$DEFAULT_BRANCH}/precice.f90"
fi

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build
  cmake --build build
fi

mkdir -p output

./build/FluidSolver ../precice-config.xml

close_log
