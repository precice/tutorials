---
title: Partitioned Pipe Multiscale
permalink: tutorials-partitioned-pipe-multiscale.html
keywords: OpenFOAM, python
summary: The 1D-3D Partitioned Pipe is a simple geometric multiscale case, that consists of flow between two pipes with heterogeneous dimensionality. The flow is incompressible and laminar. This tutorial contains an OpenFOAM case for the 3D participant and Python solver for the 1D participant.  
---

{% note %}
Get the [case files of this tutorial](https://github.com/ezonta/tutorials/tree/GeoMultiScaleTutorial/1D-3D-partitioned-pipe). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Prerequisites

- preCICE
- Python bindings of preCICE
- OpenFOAM together with the OpenFOAM adapter of preCICE

## Setup

We exchange velocity data from the 1D to the 3D participant and for the pressure data vice versa. The config looks as follows:

![Config Visualization](images/tutorials-partitioned-pipe-multiscale-config.png)

## How to run

In two different terminals execute

```bash
cd fluid1d-python && ./run.sh
```

```bash
cd fluid3d-openfoam && ./run.sh
```

## Results

You should be able to see an established laminar profile at the inlet of the 3D participant.

![Expected Result](images/tutorials-partitioned-pipe-multiscale-profile.png)
