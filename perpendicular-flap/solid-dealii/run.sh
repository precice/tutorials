#!/bin/bash
cd ${0%/*} || exit 1    		    # Run from this directory

# Solid participant

# Run this script in one terminal and the execute the Fluid participant 
# in another terminal.
# These scripts present how the two participants would be started manually.

# 1 for true, 0 for false
nonlinear=0
if [ "$1" = "-nonlinear" ]; then
            nonlinear=1
fi

linear=0
if [ "$1" = "-linear" ]; then
            linear=1
fi

if [ $linear -eq 1 ]; then
        ./linear_elasticity linear_elasticity.prm
elif [ $nonlinear -eq 1 ]; then
        ./nonlinear_elasticity nonlinear_elasticity.prm
else
        echo "No solver type specified. Please specify -linear or -nonlinear as solver type"
fi
