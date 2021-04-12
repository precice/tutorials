#!/bin/sh
set -e -u

if [ "${1:-}" = "-parallel" ]; then
    mpiexec -np 4 ./build/FluidSolver ../precice-config.xml 100 0.01 100 -parallel
else
    ./build/FluidSolver ../precice-config.xml 100 0.01 100
fi
