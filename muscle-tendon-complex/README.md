---
title: Muscle-tendon complex
permalink: tutorials-muscle-tendon-complex.html
keywords: multi-coupling, OpenDiHu, skeletal muscle
summary: In this case, an skeletal muscle (biceps) and three tendons are coupled together using a fully-implicit multi-coupling scheme.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/muscle-tendon-complex). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Case Setup

In the following tutorial we model the contraction of a muscle, in particular, the biceps. The biceps is attached to the bones by three tendons (one at the bottom and two at the top). We enforce an activation in the muscle which results in its contraction. The tendons move as a result of the muscle contraction. In this case, a muscle and three tendons are coupled together using a fully-implicit multi-coupling scheme. The case setup is shown here:

![Setup](images/tutorials-muscle-tendon-complex-setup.png) 

The muscle participant (in red), is connected to three tendons. The muscle sends traction values to the tendons, which send displacement and velocity values back to the muscle. The end of each tendon which is not attached to the muscle is fixed by a dirichlet boundary condition (in reality, it would be fixed to the bones).

The muscle and tendon meshes are obtained from patient imaging. The interfaces of the tendons and the muscle do not perfectly match, which is a quite common issue due to the limitations of imaging methods and postprocessing tools. Nonetheless, preCICE coupling methods are robust and can handle meshes that do not match perfectly. 

TODO: Explain how is the muscle activated!

## Why multi-coupling?

This is a case with four participants: the muscle and each tendon. In preCICE, there are two options to [couple more than two participants](https://www.precice.org/configuration-coupling-multi.html). The first option is a composition of bi-coupling schemes, in which we must specify the exchange of data in a participant to participant manner. However, such a composition is not suited for combining multiple strong fluid-structure interactions [1]. Thus, in this case, we use the second option, fully-implicit multi-coupling. For another multi-coupling tutorial, you can refer to the [multiple perpendicular flaps tutorial](http://precice.org/tutorials-multiple-perpendicular-flaps.html).

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

OpenDiHu is used for the muscle and each tendon participants. 
The muscle solver consists of a ... TODO
The tendon solver consists of a ... TODO

## Running the Simulation 

1. Preparation: ... TODO
   - Install OpenDiHu
   - Download input files for OpenDiHu 
   - Setup `$OPENDIHU_HOME` to your `.bashrc` file
   - Compile muscle and tendon solvers

   ```bash
   cd opendihu-solver
   ./build.sh
   ```
   - Move executables to participants directory

2. Starting:

   We are going to run each solver in a different terminal. It is important that first we navigate to the simulation directory so that all solvers start in the same directory.
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
[1] H. Bungartz, F. Linder, M. Mehl, B. Uekermann. A plug-and-play coupling approach for parallel multi-field simulations. _Comput Mech_ **55**, 1119-1129 (2015). https://doi.org/10.1007/s00466-014-1113-2

{% disclaimer %}
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®  and OpenCFD®  trade marks.
{% enddisclaimer %}
