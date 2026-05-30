#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1 || true

if [ ! -f thirdparty/precice.f90 ]; then
  # Get the preCICE Fortran module. Switch the branch with ./run.sh <branch>.
  DEFAULT_BRANCH="main"
  echo "Fetching precice.f90 from  https://github.com/precice/fortran-module/tree/${1:-$DEFAULT_BRANCH}..."
  curl --create-dirs -o thirdparty/precice.f90 "https://raw.githubusercontent.com/precice/fortran-module/${1:-$DEFAULT_BRANCH}/precice.f90"
fi

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build
  cmake --build build
fi

./build/SolidSolver ../precice-config.xml

close_log
