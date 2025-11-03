---
title: Free flow over porous media 2D
permalink: tutorials-free-flow-over-porous-media-2d.html
keywords: DuMux, porous media
summary: Flow-flow coupling example with porous media field and free flow field.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/free-flow-over-porous-media-2d). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup

This tutorial solves a coupled system consisting of a one-phase free flow and a one-phase flow in a porous media.

A pressure gradient is applied to the free flow domain from left to right. The top edge of the free-flow is a non-permeable wall with no-slip boundary conditions. In the porous media, there is a no-flow condition across the domain boundaries (left, bottom, and right boundaries). At the interface, a no-slip condition applies. The case is stationary (solved to a steady-state solution).

The setting is illustrated in the following figure:

![Free flow over porous media setup](images/tutorials-free-flow-over-porous-media-setup.png)

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-free-flow-over-porous-media-precice-config-visualization.png)

## Available solvers

Both the participants are computed using the simulation code [DuMu<sup>x</sup>](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/).

## Solver setup

To solve the flows with the DuMux framework, the necessary DUNE modules need to be downloaded and set up. This is done by running `sh setup-dumux.sh` in the tutorial folder.

Note that if an existing installation of DUNE modules is detected in a default location, this may lead to problems in running the `setup-dumux.sh` script. The script suppresses the environment variable `DUNE_CONTROL_PATH`.

To only recompile the participants, run `sh compile-dumux-cases.sh` in the tutorial folder.

## Running the simulation

Each participant has a `run.sh` script.

To run the free-flow participant, run:

```bash
cd free-flow-dumux
./run.sh
```

To run the porous-media participant, run:

```bash
cd porous-media-dumux
./run.sh
```

Participants can be executed only in serial. Parallel execution is not supported. The case takes approximately two minutes to finish.

## Post-processing

Both participants write VTU outputs, which can be viewed using ParaView.

## Further information

The results of the pressure and the velocity fields are as follows:

![Free flow over porous media results - pressure](images/tutorials-free-flow-over-porous-media-result-pressure.png)
![Free flow over porous media results - velocity](images/tutorials-free-flow-over-porous-media-result-ux.png)

Each solver folder contains an input file (`params.input`) that will be passed to the solver executables. This is a DuMUX input file describing the simulation setting, e.g., pressure, mesh size, time stepping, etc.
