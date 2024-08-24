#!/usr/bin/env bash
set -e -u

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [$1 = "Default"]; then
    echo "Running simulation with default FEniCS implementation"
    python3 ../solver-fenics/heat.py Dirichlet
else if [$1 = "HighOrder"]; then
    echo "Running simulation with higher order FEniCS implementation"
    python3 ../solver-fenics/heatHigherOrder.py Dirichlet
else if [$1 = "SDC"]; then
    echo "Running simulation with pySDC+FEniCS implementation"
    python3 ../solver-fenics/heat_pySDC.py Dirichlet
else
    echo "Invalid argument. Please provide either 'Default', 'HighOrder', or 'SDC'"
fi

close_log