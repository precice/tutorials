#!/bin/bash
set -e -u

if [ "$1" = "1d3d" ]; then
    ln -sf precice-config-1d-3d.xml precice-config.xml
    echo "Switched to 1D–3D configuration."
elif [ "$1" = "3d1d" ]; then
    ln -sf precice-config-3d-1d.xml precice-config.xml
    echo "Switched to 3D–1D configuration."
elif [ "$1" = "1d1d" ]; then
    ln -sf precice-config-1d-1d.xml precice-config.xml
    echo "Switched to 1D–1D configuration."
elif [ "$1" = "3d3d" ]; then
    ln -sf precice-config-3d-3d.xml precice-config.xml
    echo "Switched to 3D–3D configuration."
else
    echo "Usage: ./set-case.sh [1d3d|3d1d|1d1d|3d3d]"
fi

