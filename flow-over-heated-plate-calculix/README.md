---
title: Flow over heated plate with two meshes
permalink: tutorials-flow-over-heated-plate-two-meshes.html
keywords: tutorial, CHT, conjugate-heat transfer, OpenFOAM, CalculiX
summary: This tutorial describes how to run a conjugate heat transfer coupled simulation using preCICE and CalculiX as solid solver, which requires the usage of two meshes instead of one.
---


{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

The scenario is exactly the same as the one described in the [flow over heated plate tutorial](http://precice.org/tutorials-flow-over-heated-plate.html). However, this tutorial is specialized for the case when heat fluxes and temperatures live on different meshes. This is the case with CalculiX: heat fluxes are defined on face centers, while temperatures are read at nodes. This requires updating the `precice-config.xml` file to take this into account. On the fluid side, a unique mesh can still be used.

## Available solvers

By default, the fluid participant reads heat-flux values and the solid participant reads temperature values for the coupled simulation. The following participants are available:

Fluid participant:

* OpenFOAM (buoyantPimpleFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://www.precice.org/adapter-openfoam-overview.html).

Solid participant:

* CalculiX. For more information, have a look at the [CalculiX adapter documentation](http://precice.org/adapter-calculix-overview.html).

## Running the Simulation

Open two separate terminals and start the desired fluid and solid participant by calling the respective run script `run.sh` located in the participant directory. For example:

```bash
cd fluid-openfoam
./run.sh
```

and

```bash
cd solid-calculix
./run.sh
```

## Post-processing

On the OpenFOAM side, you can open the `.foam` file with ParaView, or create VTK files with `foamToVTK`. CalculiX outputs `.frd` files which can be opened with `cgx` or converted into VTK files using the converter available in the adapter repository.

CalculiX produces 1000 result files, which one can then synchronize with OpenFOAM using the "Temporal Shift Scale" filter, using a scale of 0.01. For efficient visualization, use the same scale for temperature on both outputs; we recommend from 300K to 310K.

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
