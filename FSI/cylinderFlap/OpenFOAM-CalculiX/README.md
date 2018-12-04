# Tutorial for a coupled simulation with OpenFOAM and CalculiX

# Setup

The setup for this tutorial is based on the the cylinder with a flap benchmarking case from Turek and Hron. [Include cite here]. 

In the [precice/openfoam-adapter](https://github.com/precice/openfoam-adapter). Please refer to [this wiki page](https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate) of the openfoam-adapter for details and references regarding the experimental setup.

## OpenFOAM

See [Download v5.0 | Ubuntu](https://openfoam.org/download/5-0-ubuntu/). Don't forget to also update your `~/.bashrc`! See [Download v5.0 | Ubuntu -> User Configuration](https://openfoam.org/download/5-0-ubuntu/).

## Boundary conditions
The flow in the domain is determined by a fixedValue at the inlet. From time [0-2] seconds, the flow speed is linearly increased to the a fixed value of 2m/s uniform inflow, leading to a Reynolds number of 200, where the cylinder diameter is used as the length scale. 

Note that the inflow conditions is somewhat different from the parabolic inflow condition specified at the original benchmarking case. The inflow is easy to recreate, using for example the third-party package 'swak4foam' with 'groovyBC'.


## preCICE + OpenFOAM adapter

**preCICE:** See [preCICE wiki](https://github.com/precice/precice/wiki/Building). If you have problems compiling, see the "Troubleshooting" section below.
**OpenFOAM adapter:** See [OpenFOAM adapter wiki](https://github.com/precice/openfoam-adapter/wiki/Building). If you have problems compiling, see the "Troubleshooting" section below.

To make sure that everything is working properly, you should run the following tutorial case: https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate

## Calculix

ccx 2.13 is used in the creation of this tutorial. 
Calculix sourcecode can be found at [calculix.de](http://www.calculix.de/)

## Get the OpenFOAM Fluid case

Copy the Folder `Fluid` and all its contents from https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate/buoyantPimpleFoam-laplacianFoam to this folder.

# Running

To start the coupled simulation, run the command 
'Allrun'
and to clean the case, 
'Allclean'.

Alternatively, you can start the separate solvers by running './runFluid' and './runSolid' in separate terminals. 