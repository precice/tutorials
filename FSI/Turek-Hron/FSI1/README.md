# Fluid-Structure Interaction: Turek-Hron FSI1 benchmark simulation using OpenFOAM, CalculiX and preCICE

​

This is a part of our undergraduate project titled 'Aeroelastic modelling of insect flight' carried out in the School of Mechanical Sciences, Indian Institute of Technology Goa.

​

* Authors: Nithin Adidela, Revanth Sharma, Y. Sudhakar

* email: {nithin.adiela.16003,revanth.sharma.16003,sudhakar}@iitgoa.ac.in

​

​

This case was inspired by the contribution from Derek Risseeuw (TU Delft). 



Fluid part of the simulation is performed on OpenFOAM v7. Solid part of the simulation is performed on CalculiX 2.15. Coupling of the simulation is achieved by preCICE-1.6.1

​

The  simulation is tested in serial mode in a machine running Ubuntu 18.04. The mesh used is coarse and the user is free to make the mesh finer to achieve a better agreement with the standard case after properly understanding the blockMeshdict file.

​

​

## Installing Groovy BC

​

Run the following commands in your machine, preferably in the **$FOAM_RUN** directory to download, install and use groovyBC. It takes about 30 minutes to install

​

> $ hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam -u develop\

> $cd swak4Foam\

> $./maintainanceScripts/compileRequirements.sh\

> $./Allwmake

​

After installing, edit your bash script or give an echo command to the terminal to let it discover the library. It should look something like this.

​

>  $ export SWAK4FOAM_SRC="/home/nithin/OpenFOAM/nithin-7/run/swak4Foam/Libraries:${SWAK4FOAM_SRC}"

​

​

## Running the case

​

The case files are prepared for the latest versions of OpenFOAM v7 and *pimpleFoam* solver is used. In case you are using a previous OpenFOAM version, you need to adjust the solver to *pimpleDyMFoam* in the *Fluid/system/controlDict* file. The simulations can be run in serial mode using the *Allrun* script as follows

​

> $ cd FSI1\

> $ ./Allrun

​

​

​

## Check the progress of simulations

​

The following commands can be used to check the progress

​

> $ tail -f Fluid.log\

> $ tail -f Solid.log

​

## Outputs

​

Y-displacement of the watchpoint can be plotted from start point to the present time by runnning

​

> $ ./plotDisplacement.sh

​


​

## Notes

*   CGX is used to prepare the Solid participant and it should be executable by running 'cgx'. If it has a different name (e.g. 'cgx_2.13'), adapt the respective run script, or set an alias/link from 'cgx' to 'cgx_2.13'.

*  There is an [open issue](https://github.com/precice/openfoam-adapter/issues/26) that leads to additional empty result directories when running with some OpenFOAM versions, leading to inconveniences during post-processing. Please run the script 'removeObsoleteSolvers.sh' ($ ./removeObsoleteSolvers.sh) to delete the additional files.

* You may adjust the end time in the *precice-config_*.xml*, or interupt the execution earlier if you want.


