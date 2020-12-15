---
title: Turek Hron FSI3
permalink: tutorials-turek-hron-fsi3.html
keywords: OpenFOAM, deal.II, verification
summary: The Turek-Hron FSI cases are well-established numerical benchmarks and, therefore, well suited for verification of preCICE itself and the used adapters. In this tutorial, we focus on the FSI3 case, which presents the most challenging case in terms of added mass. Please note that the meshes of this case are significantly finer than for other tutorials. Running the simulation might take a few hours. We do not recommend to run this tutorials as your first preCICE tutorial.  
---

## Setup

The setup is shown schematically here:

![FSI3 setup](images/tutorials-turek-hron-fsi3-setup.png)

For more information please refer to the original publication of the benchmark [1].

## Available solvers

Fluid participant:

* OpenFOAM. For more information, have a look at the [OpenFOAM adapter documentation](adapter-openfoam-overview.html). 

{% include important.html content="For the parabolic inflow profile, this tutorial requires groovyBC. groovyBC is part of swak4Foam. You can find more explanations in the [official OpenFOAM wiki](https://openfoamwiki.net/index.php/Contrib/swak4Foam) or get it from an [inofficial GitHub mirror](https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-swak4Foam-dev.git). Please follow the building instructions there." %}

Solid participant:

* deal.II. For more information, have a look at the [deal.II adapter documentation](adapter-dealii-overview.html). This tutorial requires the nonlinear solid solver. Please copy the nonlinear solver executable to the `solid-dealii` folder or make it discoverable at runtime and update the `solid-dealii/run.sh` script.


## Running the Simulation

Open two separate terminals and start each participant by calling the respective run script.

```
cd fluid-openfoam
./run.sh
```
and
```
cd solid-dealii
./run.sh
```

You can also run OpenFOAM in parallel by `./run.sh -parallel`. The default setting here uses 25 MPI ranks. You can change this setting in `fluid-openfoam/system/decomposeParDict`.
Moreover, the name of your solver might differ, depending on your OpenFOAM version. Have a look in the `fluid-openfoam/system/controlDict` file and set the appropriate solver name.

You may adjust the end time in the `precice-config.xml`, or interupt the execution earlier if you want.

In the first few timesteps, many coupling iterations are required for convergence. Don't lose hope, things get better quickly. 


## Postprocessing
   
You can visualize the results of the coupled simulation using e.g. ParaView. Fluid results are in the OpenFOAM format and you may load the `Fluid.foam` file. Solid results are in VTK format. If you want to visualize both domains with ParaView, keep in mind that the deal.II solver writes results every few timesteps, while the OpenFOAM solver writes in reference to simulated time. For this reason, make sure that you use compatible write intervals. You may also need to convert the OpenFOAM results to VTK (with the command `foamToVTK`). 

There is an [known issue](https://github.com/precice/openfoam-adapter/issues/26) that leads to additional "empty" result directories when running with some OpenFOAM versions, leading to inconveniences during post-processing. Please run the script `removeObsoleteSolvers.sh` to delete the additional files before importing your results in ParaView.

Moreover, we defined a watchpoint at the flap tip and can observe the perpendicular displacement with the script `plotDisplacement.sh`. 

TODO: update picture according to actual used mesh.

![FSI3 watchpoint](https://user-images.githubusercontent.com/33414590/58789906-882f7400-85ef-11e9-968a-082b33493f34.png)


## References

[1]  S. Turek, J. Hron, M. Madlik, M. Razzaq, H. Wobker, and J. Acker. Numerical simulation and benchmarking of a monolithic multigrid solver for fluid-structure interaction problems with application to hemodynamics. In H.-J. Bungartz, M. Mehl, and M. Schäfer, editors, Fluid Structure Interaction II: Modelling, Simulation, Optimization, page 432. Springer Berlin Heidelberg, 2010.

{% include disclaimer.html content="This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks." %}
