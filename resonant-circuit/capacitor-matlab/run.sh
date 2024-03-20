#!/bin/sh

# Import MATLAB-bindings path
. ../matlab-bindings-path.sh

# Add bindings to MATLAB path
export MATLABPATH=$path

# Run MATLAB code without GUI
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 matlab -nodisplay -nosplash -nodesktop -r "Solver_I;exit;"

