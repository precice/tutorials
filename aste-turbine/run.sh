#!/usr/bin/env bash

set -e
set -x

# This script assumes the ASTE binaries and python scripts are in $PATH

# Download the meshes
test -f meshes.tar.gz  || wget https://gitlab.lrz.de/precice/precice2-ref-paper-setup/-/raw/main/meshes/meshes.tar.gz

# Extract the meshes
tar -xvf meshes.tar.gz

# Calculate on fine mesh
./vtk_calculator.py 0.009.vtk x+y -t "x + y"

# Decompose both meshes to two procesors
./partition_mesh.py 0.009.vtk -n 2 -o fine_mesh
./partition_mesh.py 0.01.vtk -n 2 -o coarse_mesh

# The result directory of preciceMap needs to exist beforehand
mkdir -p mapped

# Map from the finer mesh to coarser mesh
mpirun -n 2 preciceMap -v -p A --mesh fine_mesh/fine_mesh --data "x + y" &
mpirun -n 2 preciceMap -v -p B --mesh coarse_mesh/coarse_mesh --output mapped/mapped --data "x + y"

# Join the output files together to result.vtk,
# Recovery cannot be used since GlobalID's are not exist in mapped mesh
./join_mesh.py mapped/mapped -o result.vtk

# Measure the difference between the original function and the mapped values
# Save into data array called difference
./vtk_calculator.py result.vtk x+y -t difference -it "x + y" --diff