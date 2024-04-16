---
title: Channel transport
permalink: tutorials-channel-transport.html
keywords: volume coupling, chemistry, OpenFOAM, Nutils, species, transport
summary: A CFD problem is coupled to a transport (of, e.g., a chemistry species) in a uni-directional way.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/channel-transport). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

We model a two-dimensional incompressible fluid flowing through a channel with an obstacle. The fluid problem is coupled to a simple transport problem in a uni-directional fashion. The transport problem is initialized with a circular blob close to the inflow. The density of the species is denoted with `T` (like temperature). The convected and diffused blob after 23 timesteps:

![Flap setup](images/tutorials-channel-transport-physics.png)

The behavior of the blob over the full 200 timesteps looks as follows:

<video autoplay loop>
  <source src="images/tutorials-channel-transport-animation.webm" type="video/webm">
  Animation of blob over 200 timesteps.
</video>

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-channel-transport-precice-config.png
)

## Available solvers

Fluid participant:

* Nutils. For more information, have a look at the [Nutils adapter documentation](https://www.precice.org/adapter-nutils.html). This Nutils solver requires at least Nutils v7.0.

* OpenFOAM (pimpleFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://precice.org/adapter-openfoam-overview.html).

Transport participant:

* Nutils. For more information, have a look at the [Nutils adapter documentation](https://www.precice.org/adapter-nutils.html). This Nutils solver requires at least Nutils v7.0.

## Running the simulation

Open two separate terminals and start one fluid and one transport participant by calling the respective run scripts `run.sh` located in each of the participants' directory. For example:

```bash
cd fluid-nutils
./run.sh
```

and

```bash
cd transport-nutils
./run.sh
```

## Post-processing

All solvers generate `vtk` files which can be visualized using, e.g., ParaView.
