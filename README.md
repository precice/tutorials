# preCICE Tutorials

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](http://www.precice.org/).

You may find step-by-step instructions for each case in the [preCICE wiki](https://github.com/precice/precice/wiki). *More tutorials come with each adapter* and you can also find them in the wiki.

The files in this repository are organised in the form `problem_type/geometry/solvers`. The following cases are currently provided:

* `CHT`: Conjugate Heat Transfer
   * `heat_exchanger`: A shell-and-tube heat exchanger
      * `buoyantSimpleFoam-CalculiX`
* `FSI`: Fluid-Structure Interaction
   * `flap_perp`: A flap attached on the walls of a channel, perpendicular to the flow
      * `SU2-CalculiX`
      * `OpenFOAM-CalculiX`
    * `cylinderFlap`: A cylinder with a flexible flap in a channel flow. The von Karman vortices cause the flap to oscillate. 
      * `OpenFOAM-CalculiX`