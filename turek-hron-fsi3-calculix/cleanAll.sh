#!/bin/sh

cd fluid-openfoam/
./clean.sh
cd ../solid-calculix/
./clean.sh
cd ..
rm -rf precice-run
