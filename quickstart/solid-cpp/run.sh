#!/bin/bash
cd ${0%/*} || exit 1    		    # Run from this directory

# Solid participant

# Run this script in one terminal and the execute the Fluid participant 
# in another terminal.
# These scripts present how the two participants would be started manually.

EXE=./rigid_body_solver
if [ -f "$EXE" ]; then
    $EXE
else 
    echo "Unable to locate the executable $EXE. Have a look at the README for building instructions."
fi
