#!/usr/bin/env sh
set -e -u

rm -rf precice-profiling
rm -f neumann-surrogate.log precice-Neumann-convergence.log precice-Neumann-iterations.log
rm -f surrogate.npz