---
title: 1D mass-spring system
permalink: tutorials-mass-spring-1d.html
keywords: Python, 1D
summary: We solve an oscillator with two masses in a partitioned fashion. Each mass is solved by an independent process.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/mass-spring-1d). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

The setup is taken from [1].

## Available solvers

This tutorial is only available in python. You will need to have preCICE and the python bindings installed on your system.

- *Python*: An example solver using the preCICE [Python bindings](https://www.precice.org/installation-bindings-python.html). This solver also depends on the Python libraries `numpy`, which you can get from your system package manager or with `pip3 install --user <package>`.

## Running the Simulation

### Python

Open two separate terminals and start each participant by calling:

```bash
python3 mass-spring.py MassOne
```

and

```bash
python3 mass-spring.py MassTwo
```

## Post-processing

TODO: Show how to perform convergence study. Show how to plot trajectory and inspect conservation of energy.

## References

[1] V. Schüller, B. Rodenberg, B. Uekermann and H. Bungartz, A Simple Test Case for Convergence Order in Time and Energy Conservation of Black-Box Coupling Schemes, in: WCCM-APCOM2022. [URL](https://www.scipedia.com/public/Rodenberg_2022a)
