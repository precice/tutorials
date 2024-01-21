---
title: Two-scale heat conduction
permalink: tutorials-two-scale-heat-conduction.html
keywords: Macro-micro, Micro Manager, Nutils, Heat conduction
summary: We solve a two-scale heat conduction problem with a predefined micro structure of two materials. One macro simulation is coupled to several micro simulations using the Micro Manager.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/two-scale-heat-conduction). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

This tutorial solves a heat conduction problem on a 2D domain which has an underlying micro-structure. This micro-structure changes the constituent quantities necessary for solving the problem on the macro scale. This leads to a two-scale problem with one macro-scale simulation and several micro-scale simulations.

![Case setup of two-scale-heat-conduction case](images/tutorials-two-scale-heat-conduction-macro-micro-schematic.png)

At each Gauss point of the macro domain there exists a micro simulation. The macro problem is one participant, which is coupled to many micro simulations. Each micro simulation is not an individual coupling participant, instead we use a managing software which controls all the micro simulations and their coupling via preCICE. The case is chosen from the first example case in the publication

*Bastidas, Manuela & Bringedal, Carina & Pop, Iuliu Sorin (2021), A two-scale iterative scheme for a phase-field model for precipitation and dissolution in porous media. Applied Mathematics and Computation. 396. 125933. [10.1016/j.amc.2020.125933](https://doi.org/10.1016/j.amc.2020.125933)*.

## Available solvers and dependencies

* Both the macro and micro simulations can be solved using the finite element library [Nutils](https://nutils.org/install.html) v7 or the simulation framework [DuMuX](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/).
* While using Nutils, the macro simulation is written in Python, so it requires the [Python bindings of preCICE](https://precice.org/installation-bindings-python.html).
* The [Micro Manager](https://precice.org/tooling-micro-manager-installation.html) controls all micro-simulations and facilitates coupling via preCICE. Use the [develop](https://github.com/precice/micro-manager/tree/develop) branch of the Micro Manager.
* To solve the same case written based on the DuMuX framework, the complete DuMuX environment should be installed with the provided script `setup-dumux.sh`.

## Running the simulation
#### With Nutils
Run the macro problem:

```bash
cd macro-nutils
./run.sh
```

Check the Micro Manager [configuration](https://precice.org/tooling-micro-manager-configuration.html) and [running](https://precice.org/tooling-micro-manager-running.html) documentation to understand how to set it up and launch it. There is a Python script `run-micro-problems.py` in the tutorial directory to directly run the Micro Manager. This script imports the Micro Manager, and calls its `initialize()` and `solve()` methods. The Micro Manager can be run via this script in serial or parallel. Run it in serial:

```bash
cd micro-nutils
./run.sh -s
```

Run it in parallel:

```bash
cd micro-nutils
./run.sh -p <num_procs>
```

The `num_procs` needs to fit the decomposition specified in the `micro-manager-config.json` (default: two ranks).

Even though the case setup and involved physics is simple, each micro simulation is an instance of Nutils, which usually has a moderately high computation time. If the Micro Manager is run on 2 processors, the total runtime is approximately 25-30 minutes. Do not run the Micro Manager in serial, because the runtime will be several hours.
#### With DuMuX
Be sure to install the DuMuX environment before running the simulation with:
```bash
bash setup-dumux.sh
```

Run the macro problem:

```bash
cd macro-dumux
./run.sh
```
Run the micro problem:

```bash
cd micro-dumux
./run.sh 
```

Run it in parallel:

```bash
cd micro-dumux
./run.sh -p <num_procs>
```

While running in parallel, the `num_procs` needs to fit the decomposition specified in the `micro-manager-config.json` (default: serial running). With serial running, the total simulation costs less than 2 minutes. 


## Post-processing

![Results of two-scale-heat-conduction case](images/tutorials-two-scale-heat-conduction-results.png)
#### Solution from Nutils
The participant `macro-nutils` outputs `macro-*.vtk` files which can be viewed in ParaView to see the macro concentration field. The Micro Manager uses the [export functionality](https://precice.org/configuration-export.html#enabling-exporters) of preCICE to output micro simulation data and [adaptivity related data](https://precice.org/tooling-micro-manager-configuration.html#adding-adaptivity-in-the-precice-xml-configuration) to VTU files which can be viewed in ParaView. To view the data on each micro simulation, create a Glyph on the Micro Manager VTU data. In the figure above, micro-scale porosity is shown. For a lower concentration value, the porosity increases (in the lower left corner).

![Evolving micro simulations](images/tutorials-two-scale-heat-conduction-evolving-micro-simulations.png)

The micro simulations themselves have a circular micro structure which is resolved in every time step. To output VTK files for each micro simulation, uncomment the `output()` function in the file `micro-nutils/micro.py`. The figure above shows the changing phase field used to represent the circular micro structure and the diffuse interface width.
#### Solution from DuMuX
Similar to the output data files from simulation with Nutils, the `VTU` files from macro and micro solvers could be rendered and inspected with ParaView with the mentioned method. 

The macro fields lie inside the folder `macro-dumux/build/appl`, and the micro solutions inside `micro-dumux/build/appl/precice-exports`.
