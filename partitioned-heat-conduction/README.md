---
title: Partitioned heat conduction
permalink: tutorials-partitioned-heat-conduction.html
keywords: FEniCS, Nutils, Heat conduction
summary:
---

## Setup

In this tutorial we solve a partitioned heat equation. For information on the non-partitioned case, please refer to [1, p.37ff]. The corresponding FEniCS code can be found on [github:fenics-tutorial](https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft03_heat.py).

In this tutorial the computational domain is partitioned and coupled via preCICE. The coupling roughly follows the approach described in [2].

## Available solvers and dependencies

You can either couple a solver with itself or different solvers with each other. In any case you will need to have preCICE and the python bindings installed on your system.

* `fenics`, requires you to install [FEniCS](https://fenicsproject.org/download/) and the [FEniCS-adapter](https://github.com/precice/fenics-adapter). 

* :construction: This case is still under construction. See https://github.com/precice/tutorials/issues/152. :construction: `nutils`, requires you to install [Nutils](http://www.nutils.org/en/latest/).

## Running the simulation

For choosing whether you want to run the Dirichlet-kind and a Neumann-kind participant, please provide the following commandline input:

* `-d` flag will enforce Dirichlet boundary conditions on the coupling interface.
* `-n` flag will enforce Neumann boundary conditions on the coupling interface.

For running the case, open two terminals and:

```
cd fenics
python3 heat.py -d
```

and

```
cd fenics
python3 heat.py -n
```

If you want to use Nutils for one or both sides of the setup, just `cd nutils`. The FEniCS case also supports parallel runs. Here, simply execute

```
mpirun -n <N_PROC> python3 heat.py -d
```

## Visualization

Output is written into the folder `fenics/out`. You can visualize the content with paraview by opening the `*.pvd` files. The files `HeatDirichlet` and `HeatNeumann` correspond to the numerical solution of the Dirichlet, respectively Neumann, problem, while the files with the prefix `ref` correspond to the analytical reference solution.

![Animation of the partitioned heat equation](HT_FEniCS_movie.gif)

## References

[1] Langtangen, Hans Petter, and Anders Logg. "Solving PDEs in Minutes-The FEniCS Tutorial Volume I." (2016). [pdf](https://fenicsproject.org/pub/tutorial/pdf/fenics-tutorial-vol1.pdf)
[2] Monge, Azahar, and Philipp Birken. "Convergence Analysis of the Dirichlet-Neumann Iteration for Finite Element Discretizations." Proceedings in Applied Mathematics and Mechanics (2016).
