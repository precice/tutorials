# ---
title: ASTE (Artificial Solver Testing Environment) Turbine Tutorial
permalink: tutorials-aste-turbine.html
keywords: ASTE, Testing, Turbine
summary: This tutorial is an example case for ASTE, showcasing basic features and usage of ASTE.
---

{% include note.html content="Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/turek-hron-fsi3). Read how in the [tutorials introduction](https://precice.org/tutorials.html)." %}

## Setup

This case is a simple usage scenario for. Some features offered by this case:

* Mapping from finer grid to coarse grid.
* Usage of ASTE partitioner.
* Usage of ASTE mesh joiner.
* Usage of ASTE Calculator for calculating a function on given mesh.
* Usage of ASTE.


## Running the ASTE

Run the `run.sh` script. It do the followings:

* Downloads the meshes from preCICE repository.
* Calculates function `x + y` on a finer grid.
* Partitions the coarse and fine mesh into 2 domains.
* Map the data from fine grid to coarse grid.
* Join the result meshes into a final mesh.
* Using calculator calculates difference between mapped and original function on coarse grid.

Important note :  The `run.sh` script assumes the ASTE binaries and python scripts are in $PATH
