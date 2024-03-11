---
title: Muscle-tendon complex
permalink: tutorials-muscle-tendon-complex.html
keywords: multi-coupling, OpenDiHu, skeletal muscle
summary: In this case, a skeletal muscle (biceps) and three tendons are coupled together using a fully-implicit multi-coupling scheme.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/muscle-tendon-complex). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Case Setup

In the following tutorial, we model the contraction of a muscle (the biceps). The biceps is attached to the bones by three tendons (one at the bottom and two at the top). We enforce an activation in the muscle which results in its contraction. The tendons move as a result of the muscle contraction. In this case, a muscle and three tendons are coupled together using a fully-implicit multi-coupling scheme. The case setup is shown in the following figure:

![Setup](images/tutorials-muscle-tendon-complex-setup.png)

The muscle participant (in red) is connected to three tendons. The muscle sends traction values to the tendons, which send displacement and velocity values back to the muscle. The end of each tendon which is not attached to the muscle is fixed by a dirichlet boundary condition (in reality, it would be fixed to the bones).

The muscle and tendon meshes are obtained from patient imaging. The interfaces of the tendons and the muscle do not perfectly match, which is a quite common issue due to the limitations of imaging methods and postprocessing tools. Nonetheless, preCICE coupling methods are robust and can handle meshes that do not perfectly match.

TODO: Explain how is the muscle activated!

## Why multi-coupling?

This is a case with four participants: the muscle and each tendon. In preCICE, there are two options to [couple more than two participants](https://www.precice.org/configuration-coupling-multi.html). The first option is a composition of bi-coupling schemes, in which we must specify the exchange of data in a participant-to-participant manner, limited to primarily explicit coupling schemes. However, such a composition is not suited for combining multiple strong interactions [1]. Thus, in this case, we use the second option, fully-implicit multi-coupling. For another multi-coupling tutorial, you can refer to the [multiple perpendicular flaps tutorial](http://precice.org/tutorials-multiple-perpendicular-flaps.html).

We can set this in our `precice-config.xml`:

```xml
<coupling-scheme:multi>
   <participant name="Muscle" control="yes"/>
   <participant name="Tendon-Bottom"/>
   <participant name="Tendon-Top-A"/>
   <participant name="Tendon-Top-B"/>
```

The participant that has the control is the one that it is connected to all other participants. This is why we have chosen the muscle participant for this task.

## About the solvers

We use solvers based on [OpenDiHu](https://github.com/opendihu/opendihu) for all participants.

**The muscle solver** consists of a multi-physcis multi-scale solver itself. It combines two OpenDiHu solvers in one: the *FastMonodomainSolver* and the *MuscleContractionSolver*. The two solvers are coupled using the OpenDiHu coupling tool for weak coupling.

- The [FastMonodomainSolver](https://opendihu.readthedocs.io/en/latest/settings/fast_monodomain_solver.html) models the electrochemical processes that take place in the muscle fibers, i.e, how an electrical signal propagates from the center to the extremes of the muscle fibers. The electrical signal triggers chemical reactions which lead to the contraction of sarcomeres, the smallest contraction unit in the muscle. The solver solves the so called "monodomain equation" independently for each fiber. The equation has a reaction term (small time scale) and a diffusion term (large time scale) and is solved using Strang splitting. The sarcomeres, i.e., the reaction term, are modelled using a variant of the Shorten model, specified by the CellML file `opendihu/examples/electrophysiology/input/2020_06_03_hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007.cellm`.

- The [MuscleContractionSolver](https://opendihu.readthedocs.io/en/latest/settings/muscle_contraction_solver.html) models the mechanics of the muscle. It consists of a dynamic FEM solver that models an hyperelastic active material. The active component is calculated from the active paramter $\gamma$, which ranges from 0 (no activation) to 1 (maximum activation) and is calculated in the *FastMonodomainSolver*. The material parameters are chosen as in [Heidlauf et al.](https://link.springer.com/article/10.1007/s10237-016-0772-7)

**The tendon solver** is a dynamic FEM mechanical solver. It models an hyperelastic passive material. The material parameters are chosen as in [Carniel et al.](https://pubmed.ncbi.nlm.nih.gov/28238424/)

## Running the Simulation

1. Preparation:
   - Install OpenDiHu

      In the OpenDiHu website you can find detailed [installation instructions](https://opendihu.readthedocs.io/en/latest/user/installation.html).
      We recommend to download the code from the [GitHub repository](https://github.com/opendihu/opendihu) and to run `make release_without_tests` in the parent directory.

      {% note %}
      OpenDiHu automatically downloads dependencies and installs them in the `opendihu/dependencies/` folder. You can avoid that by setting e.g., `PRECICE_DOWNLOAD = False` in the [user-variables.scons.py](https://github.com/opendihu/opendihu/blob/develop/user-variables.scons.py) before building OpenDiHu.
      {% endnote %}

   - Download input files for OpenDiHu

      OpenDiHu requires input files hosted in [Zenodo](https://zenodo.org/records/4705982) which include CellML files (containing model equations) and mesh files. Downloading these files is necessary to simulate muscles and/or tendons with OpenDiHu. You can [directly download the necessary files](https://zenodo.org/record/4705982/files/input.tgz?download=1). Extract the files and place them in the `opendihu/examples/electrophysiology/` directory.

   - Setup `$OPENDIHU_HOME` to your `.bashrc` file

      ```bash
      export OPENDIHU_HOME=/path/to/opendihu
      ```

   - Compile muscle and tendon solvers

      ```bash
      cd solver-opendihu
      ./build.sh
      ```

2. Starting the simulation:

   We are going to run each solver in a different terminal. It is important that first we navigate to the respective case directory, so that all solvers start from directories with same parent directory (see `exchange-directory` in the `precice-config.xml`).
   To start the `Muscle` participant, run:

   ```bash
   cd muscle-opendihu
   ./run.sh
   ```

   To start the `Tendon-Bottom` participant, run:

   ```bash
   cd tendon-bottom-opendihu
   ./run.sh
   ```

   To start the `Tendon-Top-A` participant, run:

   ```bash
   cd tendon-top-A-opendihu
   ./run.sh
   ```

   Finally, to start the `Tendon-Top-B` participant, run:

   ```bash
   cd tendon-top-B-opendihu
   ./run.sh
   ```

## Postprocessing... TODO

After the simulation has finished, you can visualize your results using e.g. ParaView.

## References TODO

<!-- markdownlint-configure-file {"MD034": false } -->
[1] H. Bungartz, F. Linder, M. Mehl, B. Uekermann. A plug-and-play coupling approach for parallel multi-field simulations. *Comput Mech* **55**, 1119-1129 (2015). https://doi.org/10.1007/s00466-014-1113-2

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
