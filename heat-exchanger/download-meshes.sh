#!/bin/sh

echo "This tutorial is based on a case prepared with SimScale."
echo "Since the mesh files are several MB large, we don't store them in the Git repository."
echo "This script downloads and extracts the missing files."
echo ""

echo "Downloading the Solid case..."
wget -nv -O - https://syncandshare.lrz.de/dl/fi3c9Xt5UzBc5hJvXzsLBHXn/Solid.tar.gz | tar -xzv
mv ./Solid/* solid-calculix
rm -r ./Solid
sed -i 's/Solid/\./g' solid-calculix/solid.inp

echo "Downloading and extracting the Inner-Fluid mesh..."
wget -nv -O - https://syncandshare.lrz.de/dl/fiNsYGC1DKzgio4jS5NhsXg7/polyMesh.org.tar.gz | tar -xzv -C fluid-inner-openfoam/constant
mv fluid-inner-openfoam/constant/polyMesh.org fluid-inner-openfoam/constant/polyMesh

echo "Downloading and extracting the Outer-Fluid mesh..."
wget -nv -O - https://syncandshare.lrz.de/dl/fiEZRQ8rcVWRkoyZvANim1R1/polyMesh.org.tar.gz | tar -xzv -C fluid-outer-openfoam/constant
mv fluid-outer-openfoam/constant/polyMesh.org fluid-outer-openfoam/constant/polyMesh

echo "Completed."
