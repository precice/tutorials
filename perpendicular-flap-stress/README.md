---
title: Perpendicular flap with stresses
permalink: tutorials-perpendicular-flap-stress.html
aliases:
  - /tutorials-perpendicular-flap-stress.html
keywords: G+Smo, fluid-structure interaction, FSI, OpenFOAM
summary: This tutorial is a modified version of the “perpendicular flap” tutorial coupling stress instead of force.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/develop/perpendicular-flap-stress), as continuously rendered here, or see the [latest released version](https://github.com/precice/tutorials/tree/master/perpendicular-flap-stress) (if there is already one). Read how in the [tutorials introduction](https://precice.org/tutorials.html).
{% endnote %}

## Setup

The scenario is exactly the same as the one described in the [perpendicular flap tutorial](https://precice.org/tutorials-perpendicular-flap.html). The only difference is that we use stresses instead of forces as data sent from the fluid to the solid participant. This requires changing the mapping constraint from conservative (forces) to consistent (stresses). To avoid a "write-consistent" combination, which [cannot be used in parallel](ttps://precice.org/configuration-mapping.html#restrictions-for-parallel-participants), we exchange both meshes.

## Configuration

preCICE configuration (image generated using the [precice-config-visualizer](https://precice.org/tooling-config-visualization.html)):

![preCICE configuration visualization](images/tutorials-perpendicular-flap-stress-precice-config.png)

## Available solvers

Fluid participant:

* OpenFOAM (pimpleFoam). For more information, have a look at the [OpenFOAM adapter documentation](https://precice.org/adapter-openfoam-overview.html).

Solid participant:

* G+Smo (perpendicular-flap-vertex-gismo). This solver includes both linear and nonlinear versions of the Newmark time integrator for time-dependent structural problems. The linear version iterates using a constant stiffness matrix. The nonlinear version iterates using an updated Jacobian matrix to account for material or geometric nonlinearity. By default, the solver runs in linear mode. To switch to nonlinear mode, add `--nonlinear` as option in `run.sh`. For more information, have a look at the [G+Smo adapter documentation](https://precice.org/adapter-gismo.html).

* solids4foam (OpenFOAM). A linear-elastic, small-strain solid model, using the same flap geometry and material properties as the [perpendicular flap tutorial](https://precice.org/tutorials-perpendicular-flap.html). The `interface` patch uses the solids4foam `solidTraction` boundary condition, reading the coupled traction from the `solidTraction` field that the adapter fills. This requires solids4foam v2.2 or later and an OpenFOAM-preCICE adapter with support for reading `Stress` in solid participants. For more information, have a look at the [solids4foam documentation](https://solids4foam.github.io/).

## Running the simulation

Open two separate terminals and start the desired fluid and solid participants by calling the respective run scripts `run.sh` located in the participants' directories. For example:

```bash
cd fluid-openfoam
./run.sh
```

and

```bash
cd solid-gismo
./run.sh
```

## Post-processing

On the OpenFOAM side, you can open the `.foam` file with ParaView, or create VTK files with `foamToVTK`.

On the G+Smo side, you can open the `.pvd` file located in the `solid-gismo/output` folder using ParaView. If you prefer not to plot the simulation, simply edit the `run.sh` script and remove the  `--plot` option.

As we defined a watchpoint on the 'Solid' participant at the flap tip (see `precice-config.xml`), we can plot it with gnuplot using the script `plot-displacement.sh.` You need to specify the directory of the selected solid participant as a command line argument, so that the script can pick-up the desired watchpoint file, e.g. `plot-displacement.sh solid-gismo` or `plot-displacement.sh solid-solids4foam`. The resulting graph shows the x displacement of the flap tip. You can modify the script to plot the force instead.

![Flap watchpoint](images/tutorials-perpendicular-flap-stress-displacement-watchpoint.png)

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
