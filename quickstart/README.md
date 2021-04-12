---
title: Quickstart
permalink: quickstart.html
keywords: tutorial, quickstart
summary: "Install preCICE on Linux (e.g. via a Debian package) and couple an OpenFOAM fluid solver (using the OpenFOAM-preCICE adapter) with an example rigid body solver in C++."
layout: "page"
comments: false
search: true
sidebar: nil
topnav: topnav
toc: false
---



## Start here

1. To get a feeling what preCICE does, watch a [short presentation](https://www.youtube.com/watch?v=FCv2FNUvKA8), a [longer training session](https://www.youtube.com/watch?v=FCv2FNUvKA8), or [click through a tutorial in your browser](http://run.precice.org/).
2. Get and install preCICE. For Linux systems, this is pretty easy, for macOS and Windows still possible with a bit more effort. Just pick what suits you best on [this overview page](installation-overview.html). Facing any problems? [Ask for help](community-channels.html).
    - For example, [download](https://github.com/precice/precice/releases/latest) and install our binary package for Ubuntu 20.04 (Focal Fossa) by clicking on it or using the following commands:
    ```shell
    wget https://github.com/precice/precice/releases/download/v2.2.0/libprecice2_2.2.0_focal.deb
    sudo apt install ./libprecice2_2.2.0_focal.deb
    ```
   - If you prefer to try everything in an isolated environment, you may prefer using our [demo virtual machine](installation-vm.html).
3. We will use OpenFOAM here and in many of our tutorial cases, so [install OpenFOAM](adapter-openfoam-support.html) (most compatible version: latest ESI/OpenFOAM.com).
4. Download and install the [OpenFOAM-preCICE adapter](adapter-openfoam-get.html).
5. Couple OpenFOAM with a C++ rigid body solver. [Find the case in our tutorials](https://github.com/precice/tutorials/quickstart) and keep reading. You can either `git clone` the [tutorials repository](https://github.com/precice/tutorials), or [directly download the current state](https://github.com/precice/tutorials/archive/master.zip).

## What's next?

To become a preCICE pro:

* Get an overview of the [preCICE docs](docs.html).
* See what users talk about in the [preCICE forum](https://precice.discourse.group/).
* Run [tutorials with other coupled solvers](tutorials.html).
* Watch some [preCICE videos](https://www.youtube.com/c/preCICECoupling/).
* Meet our [community](community.html).
* Find out how to [couple your own solver](couple-your-code-overview.html).
* Tell us [your story](community-projects.html).


## About the case

This tutorial deals with a fluid-structure interaction problem. The fluid part of the simulation is computed using OpenFOAM and the rigid body motion is a rigid body model (implemented in C++) with only a single degree of freedom, namely the deflection angle of the flap in the channel. The rigid body is fixed in the origin at (0,0) and the force exerted by the fluid on the rigid body structure causes an oscillatory rotation of the body. The simulation runs for 2.5 seconds. In order to gain more control over the rigid body oscillation, a rotational spring is applied at the rigid body origin. After 1.5 seconds we increase the spring constant by a factor of 8 to stabilize the coupled problem. Feel free to modify these parameters (directly in `rigid_body_solver.cpp`) and increase the simulation time (in `precice-config.xml`).

![overview](images/quickstart-setup.png)

### Building the rigid body solver
Before starting the coupled simulation, the rigid body solver needs to be built using CMake. You can run the following commands from this directory to build the `rigid_body_solver.cpp`
```
cd solid-cpp && cmake . && make
```

### Running the coupled simulation

You may run the two simulations in two different terminals and watch their output on the screen using the `run.sh` scripts (or `run.sh -parallel`, option only available for OpenFOAM) located in each participant directory. You can cleanup the simulation using `clean.sh`.


In serial, the simulation takes roughly 30 seconds to compute.

### Visualizing the results

You can visualize the simulation results of the `Fluid` participant using ParaView (use `paraFoam` to trigger the OpenFOAM native reader or load the (empty) file `Fluid.foam` into ParaView). The rigid body doesn't generate any readable output files, but the motion can be observed in the OpenFOAM data. In addition, one could visualize the coupling meshes including the exchanged coupling data. preCICE generates the relevant files during the simulation and stores them in the directory `coupling-meshes`. In order to visualize the results, load the VTK files in ParaView and apply a `Glyph` filter. Depending on the specific ParaView version, you might additionally need to disable the `ScaleArray` option by selecting `No scale array`, since the exchanged data might be inappropriate for a scaling operation. You can further add a `Warp By Vector` filter with `Displacement` to deform the coupling data. The result would look the following way:

![result](images/quickstart-result.png)
