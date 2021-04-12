#!/bin/sh
set -e -u

if [ "${1:-}" = "-parallel" ]; then
    mpiexec -np 4 ./build/SolidSolver ../precice-config.xml 100 -parallel
else
    ./build/SolidSolver ../precice-config.xml 100
fi
