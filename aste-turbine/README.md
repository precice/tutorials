---
title: ASTE (artificial solver testing environment) turbine tutorial
permalink: tutorials-aste-turbine.html
keywords: ASTE, mapping, data mapping, mapping configuration, turbine
summary: This tutorial is an example case for ASTE, where we investigate different preCICE mappings using ASTE.
---

{% include note.html content="Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/aste-turbine). Read how in the [tutorials introduction](https://precice.org/tutorials.html)." %}

If you are completely new to ASTE have a look at our [ASTE documentation](https://precice.org/tooling-aste.html). This tutorial shows how to setup a mapping between two (artificial) meshes using preCICE and ASTE in parallel. The executed mapping can be investigated in terms of accuracy as well as runtime.

## Setup

Our example consists of a turbine blade geometry, which was triangulated using different refinement levels. The mesh files are stored in the [GitLab repository](https://gitlab.lrz.de/precice/precice2-ref-paper-setup) and correspond to the mesh files used for the mapping tests of our [version 2 reference paper](https://doi.org/10.12688/openreseurope.14445.1). The mesh files are automatically downloaded when the `run.sh` script is executed. In our exemplariy setup, we map the mesh `0.006.vtk` (left side of the figure) to the mesh `0.01.vtk` (right side of the figure).

![Turbine setup](images/tutorials-aste-setup.png)

## Running the tutorial

All necessary steps in order to run the mapping setup are summarized in the `run.sh` script. Have a look at the comments in the run script in order to understand what is happening. In particular, the script executes the following steps:

1. Downloads the meshe files from the preCICE reference paper repository and stores all of them in the `meshes` directory. Note that this step is only executed when running the script for the first time. For our example setup, we use only the two meshes mentioned above, but the tutorial can easily be modified in order to employ a different mesh.
2. Generate input data for our mapping problem. As described in the [ASTE documentation](https://precice.org/tooling-aste.html#precice-aste-evaluate), we use the python script called `precice-aste-evaluate` in order to evaluate and store a test function on our mesh. Here, we select arbitrarily Franke's function, evaluate it on the input mesh `0.006.vtk` and store the data in a mesh called `input_mesh.vtu`. There are other pre-defined test functions available in the `precice-aste-evaluate` script. Use `precice-aste-evaluate --list-functions` to get a complete list of available functions and their definition.
3. In our case, we want to execute the mapping in parallel. In order to run ASTE in parallel, we need to partition our mesh files so that each rank receives its own mesh file. The python script `precice-aste-partition`, partitions the given mesh in the specified number of pieces (here two for the input mesh and two for the output mesh). Therefore, we execute ASTE with four ranks per participant.
4. Executing the actual mapping using `precice-aste-run`, which is the ASTE core module interfacing with preCICE. Here, we map the data from our fine input mesh to the coarse output. The executable stores the mesh and data per rank in the files `mapped/mapped...` with the data called `InterpolatedData`.
5. In order to join the scattered mesh files, we use the python script `precice-aste-join`. As a result, we get one large mesh file called `result.vtu`.
6. As a last step we investigate the accuracy of our applied mapping configuration, using `precice-aste-evaluate` again. This time, we use the `--diff` flag in order to compute the error between our test function and the mapped data. We store the difference data (`Error`) on the result mesh (`result.vtu`) as well, which allows us to visualize the error distribution, e.g., using `ParaView`. `precice-aste-evaluate` also prints several global error measures to the console, which looks as follows:

```bash
---[ASTE-Evaluate] INFO : Vertex count 3458
---[ASTE-Evaluate] INFO : Relative l2 error 0.0020057221573822875
---[ASTE-Evaluate] INFO : Maximum absolute error per vertex 0.010019551121144277
---[ASTE-Evaluate] INFO : Maximum signed error per vertex 0.007814315048550458
---[ASTE-Evaluate] INFO : Minimum absolute error per vertex 0.0
---[ASTE-Evaluate] INFO : Minimum signed error per vertex -0.010019551121144277
---[ASTE-Evaluate] INFO : Median absolute error per vertex 0.0009199747224788446
---[ASTE-Evaluate] INFO : 99th percentile of absolute error per vertex 0.005908856215884218
---[ASTE-Evaluate] INFO : 95th percentile of absolute error per vertex 0.004282379936870699
---[ASTE-Evaluate] INFO : 90th percentile of absolute error per vertex 0.003504421261280574
```

The information above is additionally stored in a JSON file called `result.stats.json`. The json file can be used in order to further process and store the resulting data.
