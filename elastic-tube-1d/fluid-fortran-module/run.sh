#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

git submodule update --init

if [ ! -d build ]; then
  mkdir build
  cmake -S . -B build -DPRECICE_FORTRAN_MODULE=fortran-module/precice.f90
  cmake --build build
fi

mkdir -p output

./build/FluidSolver ../precice-config.xml

close_log
