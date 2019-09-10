# Tutorial for a partitioned heat equation using FEniCS

This tutorial is described in the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-solving-the-heat-equation-in-a-partitioned-fashion-using-FEniCS).

## Waveform relaxation (experimental)

This branch of the tutorial uses the fenics-adapter waveform bindings. If you want to run this tutorial, you have to install the version of the FEniCS-adapter that can be found on [this branch](https://github.com/precice/fenics-adapter/tree/WaveformBindingsDraft). Further information on waveform relaxation in the FEniCS-adapter can be found in issue https://github.com/precice/fenics-adapter/issues/16 and milestone https://github.com/precice/fenics-adapter/milestone/2.

### Setup description and how to run the examples

For all setups linear interpolation in time is used for each (micro-)timestep. The commands for running the example cases are given.

**WR11:** this case uses no subcycling such that `dt_neumann = dt_dirichlet = T_window` (dt for each solver equals window size).

* ` python3 heat.py -d`
* ` python3 heat.py -n`

**WR12:** this case uses waveform relaxation such that `2*dt_neumann = dt_dirichlet = T_window`.

* ` python3 heat.py -d -wr 1 2`
* ` python3 heat.py -n -wr 1 2`

**WR22:** this case uses waveform relaxation such that `2*dt_neumann = 2*dt_dirichlet = T_window`.

* ` python3 heat.py -d -wr 2 1`
* ` python3 heat.py -n -wr 2 1`

### Plain Subcycling

For comparison, you can also run the case where plain subcycling is used (i.e. no interpolation between samples from (micro-)timesteps, only the sample at the end of the window is exchanged).

* Make sure to install `fenics-adapter` from master (i.e. without waveform relaxation)
* don't run `heat.py` but `heat_subcycling.py` from this folder.

### running the experiments

In the folder `experiments` you can find experiments with different timestep size and waveform configuration. If you traverse down the folder structure, you will find preCICE and adapter config files for running the experiments. All the output will be generated in that folder.

Example for running:

`cd experiments/WR12/dT0.1`
`python3 ../../../heat.py -d -wr 1 2 -dT 0.1`
`python3 ../../../heat.py -n -wr 1 2 -dT 0.1`
