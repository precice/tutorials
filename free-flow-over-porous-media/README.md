---
title: Free flow over porous medium 2D
permalink: tutorials-free-flow-over-porous-medium-2d.html
keywords: DuMux, porous medium
summary: Flow-flow coupling example with porous medium field and free flow field.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/free-flow-over-porous-media-2d). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup

This tutorial solves a simple coupled system consisting of a one-phase free flow and a one-phase flow in a porous medium.

A pressure gradient is applied to the free flow domain from left to right, while at the top of the free-flow we have a non-permeable wall with no-slip boundary conditions. In the porous medium, we assume no-flow across the domain boundaries (left, bottom and right boundary). At the interface we assume a no-slip condition.

 <!-- TODO: Add images of the setting-->

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/precice-config-visualization.png)

## Available solvers

Both the flow in free flow and porous medium can be solved using the simulation framework [DuMu<sup>x</sup>](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/).

## Solver setup

To solve the flows with the DuMux framework, the necessary DUNE modules need to be downloaded and set up. This is done by running `sh setup-dumux.sh` in the tutorial folder.

Note that if an existing installation of DUNE modules is detected in a default location, this may lead to problems in running the `setup-dumux.sh` script. The environment variable DUNE_CONTROL_PATH is suppressed by the script.

To recompile only the simulation, run `sh compile-dumux-cases.sh` in the tutorial folder.

## Running the simulation

You can find the corresponding `run.sh`script for running the case in the folders corresponding to the solvers you want to use.

### In serial

To run the free-flow participant, run:

```bash
cd free-flow-dumux
./run.sh
```

To run the porous-medium participant, run:

```bash
cd porous-medium-dumux
./run.sh
```

A serial simulation takes approximately 2 minutes to finish.

### In parallel

All participants can be run in parallel.

To run a participant in parallel, e.g. `free-flow-dumux`, run:

```bash
cd free-flow-dumux
./run.sh -p <num_procs>
```

where `<num_procs>` is the number of processes you want to use for the participant.

## Post-processing

The VTU files from both solvers could be rendered and inspected with ParaView.

## Further information

Each solver folder contains an input file (`params.input`) that will be passed to the solver executables. This is a DuMUX input file describing the simulation setting, e.g., pressure, mesh size, time stepping, etc.
