---
title: Breaking dam with flexible pillar: 2D variant
permalink: tutorials-breaking-dam-2d.html
keywords: FSI, OpenFOAM, CalculiX, nearest-projection, IQN-ILS
summary: Tutorial for an FSI simulation of a two-dimensional water column striking a flexible wall
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/breaking-dam-2d/breaking-dam-2d). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

The two-dimensional breaking dam case is a free surface problem. A large column of water comes into contact with a flexible wall, causing the wall to bend and the water to flow over and around the wall. A no slip boundary condition is applied at the bottom, the left, and the right boundary, and a
zero pressure condition at the top boundary. The test case runs for 1s with a time step size of dt = 0.005s, for a total of 200 time steps.

![domain](images/breaking-dam-2d.png)

## Available solvers

Fluid participant:

* OpenFOAM (interFoam). In case you are using a very old OpenFOAM version, you will need to adjust the solver to `interDyMFoam` in the `Fluid/system/controlDict` file. For more information, have a look at the [OpenFOAM adapter documentation](https://www.precice.org/adapter-openfoam-overview.html).

Solid participant:

* CalculiX. For more information, have a look at the [CalculiX adapter documentation](https://www.precice.org/adapter-calculix-overview.html).

## Running the simulation

You can start the simulation by running the script `./run.sh` located in each participant directory. OpenFOAM can be executed in parallel using `run.sh -parallel`. The default setting uses 4 MPI ranks.

## Post-processing

You can visualize the results using ParaView or `cgx`(for native CalculiX resul files), as usual.

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
