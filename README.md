**This branch of the tutorials repository represents the version of the code used for producing the results presented at the GAMM CSE 2018. See the corresponding presentation [here](https://mediatum.ub.tum.de/1467486). [preCICE v 1.3.0](https://github.com/precice/precice/releases/tag/v1.3.0) was used for the experiments.**

# preCICE Tutorials

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](http://www.precice.org/).

You may find step-by-step instructions for each case in the [preCICE wiki](https://github.com/precice/precice/wiki). *More tutorials come with each adapter* and you can also find them in the wiki.

The files in this repository are organised in the form `problem_type/geometry/solvers`. The following cases are currently provided:

* `CHT`: Conjugate Heat Transfer
   * `heat_exchanger`: A shell-and-tube heat exchanger
      * `buoyantSimpleFoam-CalculiX`
   * `flow-over-plate`: A channel flow over a hot plate
      * `buoyantPimpleFoam-fenics`: for running see `README.md`
* `HT`: Heat Transfer
   * `partitioned-heat`: Solving heat equation on a partitioned domain
      * `fenics-fenics`: for running see `README.md`
* `FSI`: Fluid-Structure Interaction
   * `flap_perp`: A flap attached on the walls of a channel, perpendicular to the flow
      * `SU2-CalculiX`


