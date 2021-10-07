## Setup

The solid-dune participant consists of the following files:

* Dune module containing the files to build the tutorial
* A file to run the tutorial
* A cleanup file to delete the precice and paraview output
* The mesh file containing the grid of the solid participant

The Dune module depends on the Dune core modules, the Dune-preCICE adapter and
on Dune-elastodynamics providing the necessary tools for structural simulation.

It can be found here:
https://github.com/maxfirmbach/dune-elastodynamics

The Dune module can be build with:
`<path-to-dune-common/bin/dunecontrol> --current all`

The generated executable file solid-duneis found in:
`dune-perpendicular-flap/build-cmake/src/`

and needs to be copied to the solid-dune folder.
