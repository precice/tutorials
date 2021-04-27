#!/bin/sh
set -e -u

# File check solution from: https://stackoverflow.com/questions/91368/checking-from-shell-script-if-a-directory-contains-files

# Plot diameter from fluid-cpp
if [ -n "$(ls -A ./fluid-cpp/output/*.vtk 2>/dev/null)" ]; then
    python3 plot-vtk.py diameter fluid-cpp/output/out_fluid_ &
else
    echo "No results to plot from fluid-cpp."
fi

# Plot diameter from fluid-python
if [ -n "$(ls -A ./fluid-python/output/*.vtk 2>/dev/null)" ]; then
    python3 plot-vtk.py diameter fluid-python/output/out_fluid_ &
else
    echo "No results to plot from fluid-python."
fi
