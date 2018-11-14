# Tutorial for a partitioned heat equation, using FEniCS.

## Background

This tutorial is used to solve a partitioned heat equation. For information on the non-partitioned case, please refer to [1, p.37ff]. The corresponding FEniCS code can be found on https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft03_heat.py.

In this tutorial the computational domain is partitioned and coupled via preCICE. The coupling roughly follows the approach described in [2].

## Dependencies & Running

For running this tutorial, you have to install

* FEniCS, see https://fenicsproject.org/download/
* preCICE, see https://github.com/precice/precice/wiki
* preCICE python bindings, see https://github.com/precice/precice/wiki/Non%E2%80%93standard-APIs

As soon as all dependencies are installed, you have to open two shells and run the following commands:

Run 
```
python heat.py precice-config.xml -d
```
in the first shell to start the solver that computes the Dirichlet problem (left part of domain) and run
```
python heat.py precice-config.xml -n
```
to start the solver that computes the Neumann problem (right part of domain).

## Visualization

Output is written into the folder `out`. You can visualize the content with paraview by opening the `*.pvd` files. The files `HeatDirichlet` and `HeatNeumann` correspond to the numerical solution of the Dirichlet, respectively Neumann, problem, while the files with the prefix `ref` correspond to the analytical reference solution.

![](movie.gif)

## References

[1] Langtangen, Hans Petter, and Anders Logg. "Solving PDEs in Minutes-The FEniCS Tutorial Volume I." (2016).  
[2] Monge, Azahar, and Philipp Birken. "Convergence Analysis of the Dirichlet-Neumann Iteration for Finite Element Discretizations." Proceedings in Applied Mathematics and Mechanics (2016).
