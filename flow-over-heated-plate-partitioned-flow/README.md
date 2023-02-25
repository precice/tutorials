---
title: Flow over heated plate with partitioned flow
permalink: tutorials-flow-over-heated-plate-partitioned-flow.html
keywords: tutorial, CHT, conjugate-heat transfer, OpenFOAM, FEniCS, Nutils, FF, flow partitioning
summary: This tutorial describes how to run a conjugate heat transfer coupled simulation using preCICE and any fluid-solid solver combination of our <a href="adapters-overview.html">officially provided adapter codes</a>.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate-partitioned-flow). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

The setup for this tutorial is similar to the [flow over a heated plate](https://www.precice.org/tutorials-flow-over-heated-plate.html). In this case we additionally partition the OpenFOAM fluid to create a three-way coupling using CHT (conjugate heat transfer) and FF (fluid-fluid coupling).

The test case is two-dimensional and uses a serial-implicit coupling with Quasi-Newton acceleration for the fluid-fluid coupling. The CHT coupling between the solid and the fluid2 participant is changed to serial-explicit because it does not make sense numerically to have multiple serial-implicit schemes. (see also [here](https://precice.org/configuration-coupling-multi.html))

The flow partitioning is done with the fluid-fluid module of the [preCICE OpenFOAM adapter](https://www.precice.org/adapter-openfoam-overview.html). Because we use buoyantPimpleFoam we have to tell the adapter that the coupled pressure has the name `p_rgh`. The temperature is coupled at the fluid-solid interface by exchanging *Temperature* and *Heat-Flux*, at the fluid-fluid interface it is *FlowTemperature* and *FlowTemperatureGradient*.

## Available solvers

The following participants are available:

Fluid1 and Fluid2 participant:

* OpenFOAM (buoyantPimpleFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://www.precice.org/adapter-openfoam-overview.html).

Solid participant:

* OpenFOAM (laplacianFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://www.precice.org/adapter-openfoam-overview.html).

* FEniCS. For more information, have a look at the [FeniCS adapter documentation](https://www.precice.org/adapter-fenics.html).

* Nutils. For more information, have a look at the [Nutils adapter documentation](https://precice.org/adapter-nutils.html).

* Dune-Fem. For more information, have a look at the [official documentation of Dune-Fem](https://www.dune-project.org/sphinx/dune-fem/). The solver can be installed through [PyPI](https://pypi.org/project/dune-fem/). Make sure that you are in a Python virtual environment first, which you can create inside the `solid-dune` directory and load again before running (you may need to install some tools again in this environment). Please note that Dune-Fem uses just-in-time compilation: The first time you run the solver script, it will take some time.

## Running the Simulation

All listed solvers can be used in order to run the simulation. Open three separate terminals and start the OpenFOAM fluid participant and any solid participant by calling the respective run script `run.sh` located in the participant directory. For example:

```bash
cd fluid1-openfoam
./run.sh
```

and

```bash
cd fluid2-openfoam
./run.sh
```

and

```bash
cd solid-fenics
./run.sh
```

## Post-processing

Have a look at the [flow-over heated-plate](https://www.precice.org/tutorials-flow-over-heated-plate.html) tutorial for the general aspects of post-processing.

An example of the visualized expected results looks as follows:

![result](images/tutorials-flow-over-heated-plate-partitioned-flow-results.png)

Observe that the temperature at the bottom of the plate is 310K and at the inlet 300K. On the interface, the temperature is between these values. An area of higher temperature is formed above the plate, which is shifted towards the front, driven by the flow. The temperature is coupled smoothly across the fluid-fluid coupling interface.

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
