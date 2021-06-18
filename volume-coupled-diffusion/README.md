---
title: Volume coupled diffusion
permalink: tutorials-volume-coupled-diffusion.html
keywords: FEniCS, Diffusion, Volume Coupling
summary: This tutorial illustrates volume coupling via a simple example.
---

{% include note.html content="Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/volume-coupled-diffusion). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html)." %}

## Setup

This case illustrates how to implement volume coupling in a simple toy problem. Two diffusion problems are coupled via volume coupling. One domain (the source) has constant non-zero Dirichlet boundary conditions. The Other domain (the drain) has Neumann boundary conditions and a zero Dirichlet boundary condition at the right edge of the domain. The quantity u flows from the source to the drain. This is realized via volume coupling:

![Case setup of volume-coupled-diffusion case](images/tutorials-volume-coupled-diffusion-setup.png)

## Available solvers and dependencies

You can either couple a solver with itself or different solvers with each other. In any case you will need to have preCICE and the python bindings installed on your system.

* FEniCS. Install [FEniCS](https://fenicsproject.org/download/) and the [FEniCS-adapter](https://github.com/precice/fenics-adapter).

## Running the simulation

This tutorial is for FEniCS. You can find the corresponding `run.sh` script in the folder `fenics`.

For choosing whether you want to create a source or drain, please provide the following commandline input:

* `-s` flag will create a source.
* `-d` flag will create a drain.

For running the case, open two terminals run:

```bash
cd fenics
./run.sh -s
```

and

```bash
cd fenics
./run.sh -d
```
