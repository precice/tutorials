#!/usr/bin/env bash
set -e -x

# This script assumes the ASTE binaries and python scripts are in $PATH or ASTE installed on your system

# Download the meshes
test -f meshes.tar.gz  || wget https://gitlab.lrz.de/precice/precice2-ref-paper-setup/-/raw/main/meshes/meshes.tar.gz

mkdir -p meshes

# Extract the meshes
test -f meshes/0.006.vtk -a meshes/0.01.vtk || tar -xvf meshes.tar.gz --directory meshes

# Generate input data for the mapping problem using the predefined Franke's function function
precice-aste-evaluate -m meshes/0.01.vtk -f "franke3d" -d "Franke" -o input_mesh.vtu

# Decompose both meshes to two procesors
# Choose resolution 0.01 mesh as coarse mesh and partition the mesh using a uniform algorithm
precice-aste-partition -m input_mesh.vtu -n 2 -o coarse_mesh --dir coarse_mesh --algorithm uniform
# Choose resolution 0.006 mesh as coarse mesh and partition the mesh using a meshfree algorithm
precice-aste-partition -m meshes/0.006.vtk -n 2 -o fine_mesh --dir fine_mesh --algorithm meshfree

# Create results directory of precice-aste-run
mkdir -p mapped

# Map from coarse mesh to fine mesh, start two ASTE instances, one for each participant
mpirun -n 2 precice-aste-run -p A --mesh coarse_mesh/coarse_mesh --data "Franke" &
mpirun -n 2 precice-aste-run -p B --mesh fine_mesh/fine_mesh --output mapped/mapped --data "InterpolatedData"

# Join the output files together to result.vtu
precice-aste-join -m mapped/mapped -o result.vtu --recovery fine_mesh/fine_mesh_recovery.json

# Measure the difference between the original function and the mapped values
# Save into data array called Error
precice-aste-evaluate -m result.vtu -f "franke3d" -d "Error" --diffdata "InterpolatedData" --diff --stats
