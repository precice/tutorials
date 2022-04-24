#!/usr/bin/env bash
set -e -x

# This script assumes the ASTE binaries and python scripts are in $PATH or ASTE installed on your system

# Download the meshes
test -f meshes.tar.gz  || wget https://gitlab.lrz.de/precice/precice2-ref-paper-setup/-/raw/main/meshes/meshes.tar.gz

# Extract the meshes
tar -xf meshes.tar.gz

# Calculate on fine mesh
vtk_calculator.py -m 0.009.vtk -f "eggholder3d" -d "EggHolder"

# Decompose both meshes to two procesors
# Choose resolution 0.009 mesh as fine mesh
partition_mesh.py -m 0.009.vtk -n 2 -o fine_mesh --dir fine_mesh --algorithm topology
# Choose resolution 0.01 mesh as coarse mesh
partition_mesh.py -m 0.01.vtk -n 2 -o coarse_mesh --dir coarse_mesh --algorithm topology

# The result directory of preciceMap needs to exist beforehand
mkdir -p mapped

# Map from the finer mesh to coarser mesh
mpirun -n 2 preciceMap -v -p A --mesh fine_mesh/fine_mesh --data "EggHolder" &
mpirun -n 2 preciceMap -v -p B --mesh coarse_mesh/coarse_mesh --output mapped/mapped --data "InterpolatedData"

# Join the output files together to result.vtk,
# Recovery cannot be used since GlobalID's are not exist in mapped mesh
join_mesh.py -m mapped/mapped -o result.vtk --recovery coarse_mesh/coarse_mesh_recovery.json

# Measure the difference between the original function and the mapped values
# Save into data array called difference
vtk_calculator.py -m result.vtk -f "eggholder3d" -d difference --diffdata "InterpolatedData" --diff --stats
