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
* Usage of ASTE partitioner.
* Usage of ASTE mesh joiner.
* Usage of ASTE calculator for calculating a function on a given mesh.
* Usage of ASTE.

## Running ASTE

Run the `run.sh` script. It performs the following steps:

* Downloads the meshes from preCICE repository.
* Using `vtk_calculator.py` script, calculates function `eggholder` on a finer grid.
* Using `partition_mesh.py` script, partitions the coarse and fine mesh into 2 domains.
* Using `preciceMap` executable, maps the data from the fine grid to the coarse grid.
* Using `join_mesh.py`script, joins the result meshes into a final mesh.
* Using ``vtk_calculator.py` script, calculates difference between mapped and original function on coarse grid.
