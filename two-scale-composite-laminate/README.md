# composite-multiscale

Two-scale coupled simulation of a composite structure using the preCICE coupling library. One meso-scale simulation is coupled to many micro-scale simulations. Both the scales are solved using [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/).

## Setup

The meso-scale model is a 3D beam structure which is being axially loaded.



The micro-scale model is a 3D single fibre structure.

## Solvers

### Meso scale

#### Abaqus



### Micro scale

#### Abaqus

### NASMAT

## Dependencies

* Abaqus

* preCICE: See the scripts to [build_precice with MPI](build_scripts/build_precice_v3_with_Intel_MPI.sh) and [wihtout MPI](build_scripts/build_precice_v3_without_MPI.sh).

* Micro Manager: See the script to [build the Micro Manager](build_scripts/build_micro_manager.sh).

## Running the simulation

This case is designed to be run on the [Great Lakes HPC cluster](https://arc.umich.edu/greatlakes/) at the University of Michigan. In principle the setup should work on any cluster which has access to an adequate amount of Abaqus licenses.

To run the case, submit a job via the job script:

```bash
sbatch submit_job.sbat 
```

## Post-processing
