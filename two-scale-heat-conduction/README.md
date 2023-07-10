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

This tutorial solves a heat conduction problem on a 2D domain which has an underlying micro-structure. The micro-structure makes the problem two-scale with a clear scale separation.

![Case setup of two-scale-heat-conduction case](images/macro-micro-schematic.pdf)

At each Gauss point of the macro domain there exists a micro simulation. The macro problem is one participant, which is coupled to many micro simulations. Both the macro and micro problems are solved using the finite element library [Nutils](http://www.nutils.org/en/stable/).

The case is chosen from the first example from the paper: Bastidas, Manuela & Bringedal, Carina & Pop, Iuliu Sorin (2021), A two-scale iterative scheme for a phase-field model for precipitation and dissolution in porous media. Applied Mathematics and Computation. 396. 125933. 10.1016/j.amc.2020.125933.

## Available solvers and dependencies

* Both the macro and micro simulations are solved using the finite element library [Nutils](http://www.nutils.org/en/stable/).

* Nutils code is written in Python

* The Micro Manager (`micro-manager.py`) is a controlling components which handles all micro-simulations
and facilitates coupling with the macro-simulation via preCICE. The macro-problem and Micro Manager are configured via JSON files.

## Running the simulation

The macro problem can be started using the command:

```(python)
python3 macro-nutils/macro.py
```

The Micro Manager can be directly from the terminal or imported into a Python script and then called from it. Such a script is already provided: [run-micro-problems.py](https://github.com/precice/tutorials/tree/master/two-scale-heat-conduction/run-micro-problems.py) which can be run as:

```(python)
python3 run-micro-problems.py
```

The script can also be run in parallel in the following way:

```(python)
mpirun -n <num_procs> python3 run-micro-problems.py
```
