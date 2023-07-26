---
title: Muscle contraction with fibers
permalink: 
keywords: OpenDiHu, FEBio
summary: The muscle contraction with fibers is a solid mechanics problem. The muscle participant, is a FEM solver is used to compute the deformation of the muscle due to internal active forces. The fibers participant takes care of the computation of active forces by solving the monodomain equation in each fiber. This tutorial contains OpenDiHu and FEBio variants of the FEM solver and an OpenDiHu fibers solvers.  
---

**Disclaimer:** This tutorial is built on the work of Silas Natterer et al. You can check out their work at their [repository](https://github.com/silasnatterer/bfp), including a student report. 

## Setup

## Available solvers

The muscle and fibers participants are supported in:

- *OpenDiHu*: Check [OpenDiHu's documentation](https://opendihu.readthedocs.io/en/latest/) and follow the installation instructions. 
- *FEBio*: TODO

Note: the fibers participant is only supported in OpenDiHu.

Before running the OpenDiHu participants you first have to build it. In order to follow the following instructions you must have added `OPENDIHU_HOME` to your system variables as well as the aliases for `mkorn` and `sr`.

```bash
cd muscle-opendihu
mkorn && sr
```

```bash
cd fibers-opendihu
mkorn && sr
```

## Running the Simulation

### OpenDiHu

Open two separate terminals and start each participant by calling the respective run script.

```bash
cd muscle-opendihu
. run.sh
```

and

```bash
cd fibers-opendihu
. run.sh
```

### FEBio

TODO

## Post-processing

TODO

## References

TODO
