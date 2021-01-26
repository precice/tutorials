---
title: Perpendicular flap
permalink: tutorials-perpendicular-flap.html
keywords: fluid-structure interaction, FSI, OpenFOAM, FEniCS, Nutils, deal.II, Calculix, SU2,
summary: This tutorial describes how to run a fluid-structure interaction using preCICE and any fluid-solid solver combination of our [officially provided adapter codes](adapters-overview.html).
---

## Setup

In the following tutorial we model a two-dimensional fluid flowing through a channel. A solid, elastic flap is fixed to the floor of this channel. The flap oscillates due to the fluid pressure building up on its surface. The setup is shown schematically here:

![Flap setup](images/setupDrawing.png)

The simulated flow domain is 6 units long (x) and 4 units tall (y). The flap is located at the center of the bottom (x=0) and is 1 unit long (y) and 0.1 units thick (x). We set the fluid density $\{rho}_F= 1.0kg/m^{3}$, the kinematic viscosity $\{nu}_f= 1.0m^{2}/s$, the solid density $\{rho}_s= 3.0·10^{3}kg/m^{3}$, the Young’s modulus to $E= 4.0·10^{6} kg/ms^{2}$and the Poisson ratio $\{nu}_s = 0.3$. On the left boundary a constant inflow profile in x-direction of 10m/s is prescribed. The right boundary is an outflow and the top and bottom of the channel as well as the surface of the flap are no-slip walls.

## Available solvers

Fluid participant:

* OpenFOAM. For older OpenFOAM versions, the solver name will differ. If you are using OpenFOAM v1712 / 5.x or older have a look in the `fluid-openfoam/system/controlDict` file and set the appropriate solver name. The solver can run in parallel using the command line argument `run.sh -parallel`. For more information, have a look at the [OpenFOAM adapter documentation](adapter-openfoam-overview.html).

* Nutils. For more information, have a look at the [Nutils adapter documentation](adapter-nutils-overview.html).

* SU2. As opposed to the other fluid codes, SU2 is in particular specialized for compressible flows. Therefore the default simulation parameters haven been adjusted in order to pull the setup into the compressible flow regime. For more information, have a look at the [SU2 adapter documentation](adapter-su2-overview.html).

Solid participant:

* FEniCS. The structural model is currently limited to linear elasticity. For more information, have a look at the [FeniCS adapter documentation](adapter-fenics.html).

* CalculiX. In order to allow a reasonable comparison to all solid codes, the geometrically non-linear solver has been disabled and only a linear model is used by default. For more information, have a look at the [CalculiX adapter documentation](adapter-calculix-overview.html)

* deal.II. This tutorial is intended to be used with the linear solid solver, but works with the nonlinear solver as well. Please copy the solver executables to the `solid-dealii` folder or make it discoverable at runtime and update the `solid-dealii/run.sh` script. For more information, have a look at the [deal.II adapter documentation](adapter-dealii-overview.html).

## Running the Simulation

All listed solvers can be used in order to run the simulation. Open two separate terminals and start the desired fluid and solid participant by calling the respective run script `run.sh` located in the participant directory, e.g.

```
cd fluid-openfoam
./run.sh
```
and
```
cd solid-fenics
./run.sh
```
in order to use OpenFOAM and FeniCS for this test case.


## Post-processing

How to visualize the simulation results depends on the selected solvers, as usual. Most of the solver generate `vtk` files which can visualized using, e.g., ParaView.

As we defined a watchpoint on the 'Solid' participant at the flap tip (see `precice-config.xml`), we can plot it with gnuplot using the script `plotDisplacement.sh.` You need to specify the directory of the selected solid participant as a command line argument, so that the script can pick-up the desired watchpoint file, e.g. `plotDisplacement solid-fenics`. The resulting graph shows either the displacement at the flap tip. You can modify the script to also plot the force instead.


![Flap watchpoint](images/displacement_watchpoint.png)
