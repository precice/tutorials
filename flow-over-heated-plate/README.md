---
title: Flow over heated plate
permalink: tutorials-flow-over-heated-plate.html
keywords: tutorial, CHT, conjugate-heat transfer, OpenFOAM, FEniCS, Nutils
summary: This tutorial describes how to run a conjugate heat transfer coupled simulation using preCICE and any fluid-solid solver combination of our <a href="adapters-overview.html">officially provided adapter codes</a>.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

This scenario consists of one fluid and one solid participant and it is inspired by Vynnycky et al. [1]. A fluid enters a channel with temperature $$ T_\infty $$, where it comes in contact with a solid plate. The plate is heated at its bottom and has a constant temperature of $$ T_{hot} $$.

![img](images/tutorials-flow-over-heated-plate-example.png)

The test case is two-dimensional and a serial-implicit coupling with Aitken underrelaxation is used for the coupling.

The inlet velocity is $$ u_{\infty} = 0.1 m/s $$, the inlet temperature is $$ T_{\infty} = 300K $$. The fluid and solid have the same thermal conductivities $$ k_S = k_F = 100 W/m/K $$. Further material properties of the fluid are its viscosity $$ \mu = 0.0002 kg/m/s $$ and specific heat capacity $$ c_p = 5000 J/kg/K $$. The Prandtl number $$ Pr = 0.01 $$ follows from thermal conductivity, viscosity and specific heat capacity. The solid has the thermal diffusivity $$ \alpha = 1 m^2/s $$. The gravitational acceleration is $$ g = 9.81 m/s^2 $$.

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-flow-over-heated-plate-precice-config.png)

## Available solvers

By default, the fluid participant reads heat-flux values and the solid participant reads temperature values for the coupled simulation. The following participants are available:

Fluid participant:

* OpenFOAM (buoyantFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://www.precice.org/adapter-openfoam-overview.html).

* SU2. For more information, have a look at the [SU2 adapter docmentation](https://www.precice.org/adapter-su2.html).

Solid participant:

* OpenFOAM (laplacianFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://www.precice.org/adapter-openfoam-overview.html).

* FEniCS. For more information, have a look at the [FeniCS adapter documentation](https://www.precice.org/adapter-fenics.html).

* Nutils. For more information, have a look at the [Nutils adapter documentation](https://precice.org/adapter-nutils.html).

* Dune-Fem. For more information, have a look at the [official documentation of Dune-Fem](https://www.dune-project.org/sphinx/dune-fem/). The `run.sh` script installs the solver from [PyPI](https://pypi.org/project/dune-fem/) into a Python virtual environment. Please note that Dune-Fem uses just-in-time compilation: The first time you run the solver script, it will take some time.

It is also possible to use CalculiX as solid solver. In that case, two coupling meshes are needed: CalculiX read/writes temperatures on nodes, but read/writes heat-fluxes on face centers. This requires some adaptation of the `precice-config.xml` file, and [a separate tutorial](tutorials-flow-over-heated-plate-two-meshes.html) has been designed for it.

## Running the Simulation

All listed solvers can be used in order to run the simulation. Open two separate terminals and start the desired fluid and solid participant by calling the respective run script `run.sh` located in the participant directory. For example:

```bash
cd fluid-openfoam
./run.sh
```

and

```bash
cd solid-fenics
./run.sh
```

in order to use OpenFOAM and FEniCS for this test case. Feel free to try different combinations, they should all run and give approximately similar results.

## Post-processing

How to visualize the simulation results depends on the selected solvers. Most of the solvers generate VTK files which can visualized using ParaView or similar tools.
In case of OpenFOAM, you can open the `.foam` file with ParaView, or create VTK files with `foamToVTK`.

An example of the visualized expected results looks as follows:

![result](images/tutorials-flow-over-heated-plate-result-openfoam.png)

Observe that the temperature at the bottom of the plate is 310K and at the inlet 300K. On the interface, the temperature is between these values. An area of higher temperature is formed above the plate, which is shifted towards the front, driven by the flow.

You may use additional filters, such as the Calculator and the Plot Over Line, to obtain the distribution of the non-dimensional temperature $$ (T-T_{inlet})/(T_{solid}-T_{inlet}) $$:

![graph](images/tutorials-flow-over-heated-plate-graph-result.png)

### Compare results

First generate the output for each run by adding export to the participant `Solid` in `precice-config.xml`:

```xml
<participant name="Solid">
  <export:vtk directory="precice-exports" />
  <receive-mesh name="Fluid-Mesh" from="Fluid" />
  <provide-mesh name="Solid-Mesh" />
  ...
</participant>
```

After that running a case from this tutorial will export data into `solid-*/precice-exports`. To visualize and compare these results run `python3 plot-final-interface-temperature.py` (You can install the required python packages by running `pip3 install -r plot-final-interface-temperature-requirements.txt`). This will plot the dimensionless temperature `theta = (T-300)/(310-300)` (with `T` being the temperature) across the coupling interface, i.e. where the solid and the fluid meet and exchange heat. The x-axis shows the x coordinate and the y-axis the dimensionless temperature `theta` at the interface. If you want to exclude certain cases, simply comment out the corresponding lines in the script. For reference see below:

![Comparison of the results with different solvers](images/tutorials-flow-over-heated-plate-results-comparison.png)

## References

[1]  M. Vynnycky, S. Kimura, K. Kanev, and I. Pop. Forced convection heat transfer from a flat plate: the conjugate problem. International Journal of Heat and Mass Transfer, 41(1):45 – 59, 1998.

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
