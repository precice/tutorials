# preCICE Tutorials

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](http://www.precice.org/).

You may find step-by-step instructions for each case in the [preCICE wiki](https://github.com/precice/precice/wiki). *More tutorials come with each adapter* and you can also find them in the wiki.

The files in this repository are organized in the form `problem_type/geometry/solvers`. The following cases are currently provided:

* `CHT`: Conjugate Heat Transfer
   * `heat_exchanger`: A shell-and-tube heat exchanger
      * `buoyantSimpleFoam-CalculiX`
   * `flow-over-plate`: A channel flow over a hot plate
      * `buoyantPimpleFoam-fenics`: for running see `README.md`
      * `buoyantPimpleFoam-nutils`: for running see `README.md`
* `HT`: Heat Transfer
   * `partitioned-heat`: Solving heat equation on a partitioned domain
      * `fenics-fenics`: for running see `README.md`
* `FSI`: Fluid-Structure Interaction
   * `flap_perp`: A flap attached on the walls of a channel, perpendicular to the flow
      * `SU2-CalculiX`
      * `OpenFOAM-CalculiX`
      * `OpenFOAM-deal.II`: for running see the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-deal.II-and-OpenFOAM)
    * `flap_perp_2D`: 2D version of flap_perp
      * `OpenFOAM-deal.II`: for running see the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-deal.II-and-OpenFOAM)
    * `cylinderFlap`: A cylinder with a flexible flap in a channel flow. The von Karman vortices cause the flap to oscillate. 
      * `OpenFOAM-CalculiX`
      * `OpenFOAM-deal.II`: for running see the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-deal.II-and-OpenFOAM)
    * `cylinderFlap_2D`: 2D version of flap_perp
      * `OpenFOAM-deal.II`: for running see the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-deal.II-and-OpenFOAM)
    * `3D_Tube`: A 3D expanding tube with a pressure inlet boundary condition.
      * `OpenFOAM-CalculiX`
* `SSI`: Structure-Structure Interaction
   * `loaded_beam`: A partitioned elastic beam fixed at both ends. A mechanical force is applied to one partition.
      * `CalculiX-CalculiX`     
