---
title: ASTE (Artificial Solver Testing Environment) Turbine Tutorial
permalink: tutorials-aste-turbine.html
keywords: ASTE, Testing, Turbine
summary: This tutorial is an example case for ASTE, showcasing basic features and usage of ASTE.
---

{% include note.html content="Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/aste-turbine). Read how in the [tutorials introduction](https://precice.org/tutorials.html)." %}

## Setup

This case is a simple usage scenario for ASTE. Some features demonstrated by this tutorial:

* Mapping from a fine grid (`0.009.vtk`) to a coarse grid (`0.01.vtk`).
* Usage of `precice-aste-partition` executable.
* Usage of `precice-aste-join` executable.
* Usage of `precice-aste-evaluate` executable for calculating a function on a given mesh.
* Usage of `precice-aste-run` executable.

## Running ASTE

Run the `run.sh` script. It performs the following steps:

* Downloads the meshes from preCICE repository.
* Using `precice-aste-evaluate` script, calculates function `eggholder3d` on a finer grid.
* Using `precice-aste-partition` script, partitions the coarse and fine mesh into 2 domains.
* Using `precice-aste-run` executable, maps the data from the fine grid to the coarse grid.
* Using `precice-aste-join`script, joins the result meshes into a final mesh.
* Using `precice-aste-evaluate` script, calculates difference between mapped and original function on coarse grid.
