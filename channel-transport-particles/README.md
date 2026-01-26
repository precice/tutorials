---
title: Channel transport with particles
permalink: tutorials-channel-transport-particles.html
keywords: volume coupling, particles, OpenFOAM, MercuryDPM, transport, just-in-time mapping
summary: A CFD problem is coupled to a particles in a uni-directional way for particles tracing.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/channel-transport-particles). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup

We model a two-dimensional incompressible fluid flowing through a channel with an obstacle. The fluid problem is coupled to a particle participant for particle tracing. Similar to the transport problem (see the [channel-transport tutorial](tutorials-channel-transport.html)), particles are arranged in a circular blob close to the inflow.

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-channel-transport-precice-config.png
)

## Available solvers

Fluid participant:

* Nutils. For more information, have a look at the [Nutils adapter documentation](https://precice.org/adapter-nutils.html). This Nutils solver requires at least Nutils v7.0.

* OpenFOAM (pimpleFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://precice.org/adapter-openfoam-overview.html).

Particle participant:

* MercuryDPM

## Running the simulation

For the fluid solver, use Nutils for ease of installation and OpenFOAM for speed.

Open two separate terminals and start one fluid and one transport participant by calling the respective run scripts `run.sh` located in each of the participants' directory. For example:

```bash
cd fluid-nutils
./run.sh
```

and either the non-adaptive mesh transport solver

```bash
cd particles-mercurydpm
./run.sh
```

## Post-processing

All solvers generate `vtk` files which can be visualized using, e.g., ParaView.
