---
title: Growing Mesh
permalink: tutorials-growing-mesh.html
keywords: python, remeshing
summary: The growing mesh case is a showcase example of two solvers which grow their mesh at predefined points in time.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/growing-mesh). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup

The problem consists of a unit-square uniformly discretized by 768 x 768 nodes.
Running in parallel is only allowed for 1, 4, 9, or 16 ranks.
The unit square is partitioned equally among the ranks of a solver.

The mesh starts with 2 nodes in z direction and at a given frequency, 2 nodes are added to the mesh, changing only the load per rank, not the partitioning.

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-growing-mesh-precice-config.png)

## Available solvers

There are two solvers that define the same mesh:

- A who runs first
- B who runs second

## Running the Simulation

Pass the amount of ranks to the run script of the solvers.
Not passing a number, runs the simulation on a single rank.
To run both on a two rank each, use:

```bash
cd A
./run.sh 2
```

and

```bash
cd B
./run.sh 2
```
