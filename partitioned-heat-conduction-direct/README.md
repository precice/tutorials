---
title: Partitioned heat conduction (direct access setup)
permalink: tutorials-partitioned-heat-conduction-direct.html
keywords: Nutils, Heat conduction, Direct mesh access
summary: This tutorial is a modified version of the "partitioned heat conduction" tutorial showcasing direct mesh access.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/develop/partitioned-heat-conduction-direct), as continuously rendered here, or see the [latest released version](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-direct) (if there is already one). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup

This case is a modified version of the [partitioned heat conduction tutorial](tutorials-partitioned-heat-conduction.html). Main modification is that we here use the [direct mesh access feature](couple-your-code-direct-access.html) to let the solvers compute the data mapping and not preCICE.

Further minor modifications:

- We use a parallel coupling scheme instead of a serial one to prevent running into the problem where we are trying to add a zero column to the quasi-Newton matrix. For serial coupling, this happens here because one data field converges much faster than the other.

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-partitioned-heat-conduction-direct-precice-config.png)

The data mapping is computed by directly sampling the FEM function representation at the inquired locations.

## Available solvers

You can either couple a solver with itself or different solvers with each other. Available solvers:

- Nutils. For more information, have a look at the [Nutils adapter documentation](https://precice.org/adapter-nutils.html).

- G+Smo. Install the [G+Smo adapter](https://precice.org/adapter-gismo.html).

## Running the simulation

Open two terminals and run:

```bash
cd neumann-nutils
./run.sh
```

and

```bash
cd dirichlet-nutils
./run.sh
```

See the [partitioned heat conduction tutorial](tutorials-partitioned-heat-conduction.html).
