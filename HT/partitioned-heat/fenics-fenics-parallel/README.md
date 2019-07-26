# Tutorial for a partitioned heat equation using FEniCS

This tutorial is described in the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-solving-the-heat-equation-in-a-partitioned-fashion-using-FEniCS).

## Notes on parallel fenics runs

run: `mpirun -np 4 python3 NAME_OF_SCRIPT`

## Results for parallel run

If we want to run the simulation on a total of 6 processes, we can use the following example setup:

```
mpirun -np 4 python3 heat.py -d
mpirun -np 2 python3 heat.py -n
```

We obtain the exact solution on the following partitioning:

![Partitioning for 4 nodes for the Dirichlet problem and 2 nodes for the Neumann problem.](partitioning.png)
