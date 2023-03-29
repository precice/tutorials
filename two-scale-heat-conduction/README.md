# coupled-heat-conduction

<a style="text-decoration: none" href="https://github.com/precice/fenics-adapter/blob/master/LICENSE" target="_blank">
    <img src="https://img.shields.io/github/license/IshaanDesai/coupled-heat-conduction.svg" alt="GNU LGPL license">
</a>

This code solves a heat conduction problem on a 2D domain which has an underlying micro-structure. The micro-structure makes the problem two-scale with a clear scale separation.
At each Gauss point of the macro-domain there exists a micro-simulation. The macro-domain is resolved in the file `macro-heat.py`
and the micro-domain is resolved in the file `micro_sim/micro_heat_circular.py`. Both the macro and micro problems are solved using the finite element library [Nutils](http://www.nutils.org/en/stable/).

The coupling between the macro-simulation and several micro-simulations is achieved using the coupling library [preCICE](https://precice.org/) 
and a Micro Manager. The Micro Manager (`micro-manager.py`) is a controlling components which handles all micro-simulations
and facilitates coupling with the macro-simulation via preCICE. The macro-problem and Micro Manager are configured via JSON files.

The case is chosen from the first example from the paper: Bastidas, Manuela & Bringedal, Carina & Pop, Iuliu, (2021), A two-scale iterative scheme for a phase-field model for precipitation and dissolution in porous media. Applied Mathematics and Computation. 396. 125933. 10.1016/j.amc.2020.125933. 

## Dependencies

* **Nutils** can be installed through the [installation procedure](http://www.nutils.org/en/latest/intro/#installation).
* **preCICE** can be installed in [several ways](https://precice.org/installation-overview.html).
* [pyprecice]()
* [micro-manager]()

## Running two-scale coupled heat conduction problem

The coupled macro problem can be started using the command:

```(python)
python3 macro-heat.py
```

For a coupled simulation the micro problems are managed by the micro manager. The micro-manager is imported into a Python script and then called from it. In this case, the script is [run-micro-problems.py]() which can be run as:

```(python)
python3 run-micro-problems.py
```

The script can also be run in parallel in the following way:

```(python)
mpirun -n <num_procs> python3 run-micro-problems.py
```
