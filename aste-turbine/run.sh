#!/usr/bin/env bash
set -e -x

# This script assumes the ASTE binaries and python scripts are in $PATH or ASTE installed on your system

# Download the meshes
test -f meshes.tar.gz  || wget https://gitlab.lrz.de/precice/precice2-ref-paper-setup/-/raw/main/meshes/meshes.tar.gz

mkdir -p meshes

# Extract the meshes
test -f meshes/0.006.vtk -a meshes/0.01.vtk || tar -xvf meshes.tar.gz --directory meshes

# Generate input data for the mapping problem using the predefined Franke's function function
precice-aste-evaluate -m meshes/0.006.vtk -f "franke3d" -d "Franke" -o input_mesh.vtu

# Decompose both meshes to two procesors
# Choose resolution 0.006 mesh as fine mesh and partition the mesh using a uniform algorithm
precice-aste-partition -m input_mesh.vtu -n 2 -o fine_mesh --dir fine_mesh --algorithm uniform
# Choose resolution 0.01 mesh as coarse mesh and partition the mesh using a meshfree algorithm
precice-aste-partition -m meshes/0.01.vtk -n 2 -o coarse_mesh --dir coarse_mesh --algorithm meshfree

# The result directory of precice-aste-run needs to exist beforehand
mkdir -p mapped

# Map from the finer mesh to coarser mesh
mpirun -n 2 precice-aste-run -p A --mesh fine_mesh/fine_mesh --data "Franke" &
mpirun -n 2 precice-aste-run -p B --mesh coarse_mesh/coarse_mesh --output mapped/mapped --data "InterpolatedData"

# Join the output files together to result.vtu,
# Recovery cannot be used since GlobalID's are not exist in mapped mesh
precice-aste-join -m mapped/mapped -o result.vtu --recovery coarse_mesh/coarse_mesh_recovery.json

# Measure the difference between the original function and the mapped values
# Save into data array called error
precice-aste-evaluate -m result.vtu -f "franke3d" -d "Error" --diffdata "InterpolatedData" --diff --stats
