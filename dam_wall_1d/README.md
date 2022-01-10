---
title: Elastic Dam Wall 1D
permalink: tutorials-elastic-dam-wall-1d.html
keywords: FSI, OpenFOAM, CalculiX, nearest-projection, ILS
summary: Tutorial for an FSI simulation of a one-dimensional dam wall scenario
---

{% include note.html content="Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/elastic-dam-wall-1d). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html)." %}

## Setup

The 1D dam wall is a multiphase test case, where a large body of water impacts a flexible barrier/wall. The test case is modelled after 


## Available solvers

Fluid participant:

* OpenFOAM. This tutorial is known to work with OpenFOAM 7.0 with `interFoam`.

Solid participant:

* CalculiX. This tutorial is known to work with CalculiX 2.17, but it should also work with older and newer versions. For more information, have a look at the [CalculiX adapter documentation](https://www.precice.org/adapter-calculix-overview.html).

## Running the simulation

You can start the simulation by running the script `./run.sh` located in each participant directory. OpenFOAM can be executed in parallel by using an additional `run.sh -parallel` flag. The default setting uses 4 MPI ranks.

## Post-processing

You can visualize the results using paraView or `cgx`(for native CalculiX resul files), as usual. 

{% include disclaimer.html content="This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks." %}
