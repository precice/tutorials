#!/usr/bin/env bash
set -e -u

python3 -m venv --system-site-packages ../solver-fenics/.venv
. ../solver-fenics/.venv/bin/activate
pip install -r ../solver-fenics/requirements.txt

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

if [ $# -eq 0 ] 
then
    echo "Running simulation with default FEniCS implementation"
    python3 ../solver-fenics/heat.py Dirichlet
else
    case "$1" in
        -[hH]|--help)
            echo "Usage: $0 [highorder|sdc|*]"
            echo ""
            echo "      highorder: Run simulation with higher order FEniCS implementation"
            echo "      sdc      : Run simulation with pySDC+FEniCS implementation"
            echo "      *        : For every other input run the simulation with default FEniCS implementation"
            exit 0
            ;;
        highorder)
            echo "Running simulation with higher order FEniCS implementation"
            python3 ../solver-fenics/heatHigherOrder.py Dirichlet
            ;;
        sdc)
            echo "Running simulation with pySDC+FEniCS implementation"
            python3 ../solver-fenics/heat_pySDC.py Dirichlet
            ;;
        *)
            echo "Running simulation with default FEniCS implementation"
            python3 ../solver-fenics/heat.py Dirichlet
            ;;
    esac
fi

close_log
