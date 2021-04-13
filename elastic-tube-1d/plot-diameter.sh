#!/bin/sh
set -e -u

# Plot diameter from fluid-cpp
if [ -d "./fluid-cpp/output/" ]; then
    python3 plot-vtk.py diameter fluid-cpp/output/out_fluid_ &
else
    echo "No results to plot from fluid-cpp."
fi

# Plot diameter from fluid-python
if [ -d "./fluid-python/output/" ]; then
    python3 plot-vtk.py diameter fluid-python/output/out_fluid_ &
else
    echo "No results to plot from fluid-python."
fi
