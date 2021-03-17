---
title: Flow over heated plate
permalink: tutorials-flow-over-heated-plate.html
keywords: tutorial, CHT, conjugate-heat transfer, OpenFOAM, FEniCS, Nutils
summary: This tutorial describes how to run a conjugate heat transfer coupled simulation using preCICE and any fluid-solid solver combination of our [officially provided adapter codes](adapters-overview.html).
---

{% include important.html content="We have not yet ported the documentation of the preCICE tutorials from the preCICE wiki to here. Please go to the [preCICE wiki](https://github.com/precice/precice/wiki#2-getting-started---tutorials)" %}


## Setup

This scenario consists of one fluid and one solid participant and it is inspired by Vynnycky et al. [1]. A fluid enters a channel with temperature `T_\infty`, where it comes in contact with a solid plate. The plate is heated at its bottom and has a constant temperature of `T_hot`.

![img](images/tutorials-flow-over-heated-plate-example.png)

The test case is two-dimensional and a serial-implicit coupling with Aitken underrelaxation is used for the coupling.

The inlet velocity is `u_\infty = 0.1 m/s`, the inlet temperature is `T_\infty = 300K`. The fluid and solid have the same thermal conductivities `k_S = k_F = 100 W/m/K`. Further material properties of the fluid are its viscosity `mu = 0.0002 kg/m/s` and specific heat capacity `c_p = 5000 J/kg/K`. The Prandtl number `Pr = 0.01` follows from thermal conductivity, viscosity and specific heat capacity. The solid has the thermal diffusivity `\alpha = 1 m^2/s`. The gravitational acceleration is `g = 9.81 m/s^2`.

## Available solvers

By default, the fluid participant reads heat-flux values and the solid participant reads temperature values for the coupled simulation. The following participants are available:

Fluid participant:

* OpenFOAM (buoyantPimpleFoam). For more information, have a look at the [OpenFOAM adapter documentation](adapter-openfoam-overview.html).

Solid participant:

* OpenFOAM (laplacianFoam). For more information, have a look at the [OpenFOAM adapter documentation](adapter-openfoam-overview.html).

* FEniCS. For more information, have a look at the [FeniCS adapter documentation](adapter-fenics.html).

## Running the Simulation

All listed solvers can be used in order to run the simulation. Open two separate terminals and start the desired fluid and solid participant by calling the respective run script `run.sh` located in the participant directory. For example:

```
cd fluid-openfoam
./run.sh
```
and
```
cd solid-fenics
./run.sh
```
in order to use OpenFOAM and FEniCS for this test case.

## Post-processing

How to visualize the simulation results depends on the selected solvers. Most of the solvers generate `vtk` files which can visualized using, e.g., ParaView.
An example of the visualized expected results looks as follows:

![result](images/result-openfoam.png)

Observe that the temperature at the bottom of the plate is 310K and at the inlet 300K. On the interface, the temperature is between these values. An area of higher temperature is formed above the plate, which is shifted towards the front, driven by the flow.

You may use additional filters, such as the Calculator and the Plot Over Line, to obtain the distribution of the non-dimensional temperature (T-T_inlet)/(T_solid-T_inlet):

![graph](images/graph-result.png)

## References

[1]  M. Vynnycky, S. Kimura, K. Kanev, and I. Pop. Forced convection heat transfer from a flat plate: the conjugate problem. International Journal of Heat and Mass Transfer, 41(1):45 – 59, 1998.

{% include disclaimer.html content="This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks." %}
