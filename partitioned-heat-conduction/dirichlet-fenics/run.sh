#!/usr/bin/env bash
set -e -u

python3 -m venv --system-site-packages .venv
. .venv/bin/activate
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
            echo "      irk : Run simulation with higher-order implicit Runge-Kutta schemes FEniCS implementation"
            echo "      sdc : Run simulation with pySDC+FEniCS implementation"
            echo "      *   : For every other input run the simulation with default FEniCS implementation"
            exit 0
            ;;
        irk)
            echo "Running simulation with higher-order implicit Runge-Kutta schemes FEniCS implementation"
            python3 ../solver-fenics/heatHigherOrder.py Dirichlet
            ;;
        sdc)
            # install pySDC + its dependencies only if needed
            pip install git+https://github.com/Parallel-in-Time/pySDC@5.5.0
            pip install pySDC~=5.5
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
