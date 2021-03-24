---
title: Flow over heated plate steady state
permalink: tutorials-flow-over-heated-plate-steady-state.html
keywords: CHT, steady-state, Code_Aster, OpenFOAM
summary: Using a steady-state OpenFOAM solver for a CHT coupling with Code_Aster.
---

{% include important.html content="We have not yet ported the documentation of the preCICE tutorials from the preCICE wiki to here. Please go to the [preCICE wiki](https://github.com/precice/precice/wiki#2-getting-started---tutorials)" %}

## Setup

The setup for this tutorial is similar to the [flow over a heated plate](tutorials-flow-over-heated-plate.html) using OpenFOAM. In this tutorial OpenFOAM is used as the solver for the fluid domain, and Code-Aster is the solver for the solid domain. A difference here is that we are using a steady-state OpenFOAM solver for demonstration purposes, therefore the results between the two tutorials are not comparable.


## Available solvers

Fluid participant:

* OpenFOAM. We use buoyantSimpleFoam instead of the transient buoyantPimpleFoam. For more information, have a look at the [OpenFOAM adapter documentation](adapter-openfoam-overview.html).

Solid participant:

* Code_Aster. The [Code_Aster adapter documentation](adapter-code_aster.html) is oriented on this tutorial case. In particular the described configuration settings.

## Running the Simulation

Open two separate terminals and start each participant by calling the respective `run.sh` script.

## Post-processing

Firstly, enable the ParaViS view in Salome-Meca by selecting the icon in the top of the screen.

For visualizing the results of the fluid solver, go to `File -> Open ParaView File` and select the `.foam` file. If you're asked to choose a reader, please select `OpenFOAMReader` and click `Apply` to visualize the result.

For visualizing the result of the solid solver, press again `Open ParaView File` and select the `output-..rmed` group. Again, click `Apply` to visualize the result. After setting the temperature scale for both domains to 300-310 K, the following result is given for timestep 200:

![post-processing](images/tutorials-flow-over-heated-plate-steady-state-post-processing.png)

{% include disclaimer.html content="This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks." %}

{% include disclaimer.html content="This offering is not approved or endorsed by Électricité de France (EDF), producer and distributor of the Code_Aster software via www.code-aster.org, and owner of the Code_Aster trademark." %}
