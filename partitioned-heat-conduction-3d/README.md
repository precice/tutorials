---
title: 3D partitioned heat conduction
permalink: tutorials-partitioned-heat-conduction-3d.html
keywords: FEniCSx, Heat conduction
summary: We solve a simple heat equation on a 3D domain. The domain is partitioned and the coupling is established in a Dirichlet-Neumann fashion.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/develop/partitioned-heat-conduction-3d), as continuously rendered here, or see the [latest released version](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-3d) (if there is already one). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup


We solve a partitioned heat equation. For information on the two dimensional non-partitioned case, please refer to [1, p.37ff]. In this tutorial the computational domain is partitioned and coupled via preCICE. The coupling roughly follows the approach described in [2].

![Case setup of partitioned-heat-conduction case](images/tutorials-partitioned-heat-conduction-setup.png)

The computational domain can be seen in the picture above. To get a better overview of the domain, the domains of the participants are translated.
The domain of the Neumann participant is the small box, and the Dirichlet participant's domain is the rest. In the simulation, the small box is fully contained in the Dirichlet domain. This means, the coupling interface consists of five sides of the box.

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-partitioned-heat-conduction-precice-config.png)

## Available solvers and dependencies

You can either couple a solver with itself or different solvers with each other. In any case you will need to have preCICE and the python bindings installed on your system.

* FEniCSx. Install [FEniCS](https://fenicsproject.org/download/) and the [FEniCSx-adapter](https://github.com/precice/fenicsx-adapter). The code is adapted from the existing [fenics-tutorial](https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft03_heat.py) from [1].

## Running the simulation

You can find the corresponding `run.sh` script for running the case in the folders corresponding to the participant you want to use:

```bash
cd dirichlet-fenicsx
./run.sh
```

and

```bash
cd neumann-fenicsx
./run.sh
```

## Visualization

Output is written into the folders of the fenicsx solvers (`dirichlet-fenicsx/output-DIRICHLET.bp` and `dirichlet-fenicsx/output-DIRICHLET.bp`).

It is sufficient to import the folders to paraview to get the visualization of the simulation.

## References

[1] Hans Petter Langtangen and Anders Logg. "Solving PDEs in Minutes-The FEniCS Tutorial Volume I." (2016). [pdf](https://fenicsproject.org/pub/tutorial/pdf/fenics-tutorial-vol1.pdf)  