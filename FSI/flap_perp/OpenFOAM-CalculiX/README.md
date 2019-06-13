# Tutorial for an FSI simulation of an elastic flap perpendicular to a channel flow

This tutorial is described in the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-OpenFOAM-and-CalculiX).

It is known to work with OpenFOAM 4.1, 5.0, v1712, v1806 and CalculiX 2.13 with CGX, but it should also work with newer versions. Have a look on issues concerning [OpenFOAM 6 / dev](https://github.com/precice/openfoam-adapter/issues/21) and [OpenFOAM v1812](https://github.com/precice/openfoam-adapter/issues/59).

The case files are prepared for the latest versions of OpenFOAM and use the solver `pimpleFoam`. **In case you are using a previous OpenFOAM version** you need to adjust the solver to `pimpleDyMFoam` in the `Fluid/system/controlDict` file.

You may run the coupled simulation in serial using the script `Allrun` or in parallel with `Allrun -parallel` (`Allrun_parallel` is a shortcut to this). The output of each step will be redirected to log files. You can cleanup the simulation using `Allclean`.

If you prefer to run the two simulations in two different terminals and watch their output on the screen, use the (simpler) scripts `runFluid` (or `runFluid -parallel`) and `runSolid`. Please always run the script `runFluid` first.

CGX is used to prepare the Solid participant and it should be executable by running `cgx`. If it has a different name (e.g. `cgx_2.13`), adapt the respective run script, or set an alias/link from `cgx` to `cgx_2.13`.

There is an [open issue](https://github.com/precice/openfoam-adapter/issues/26) that leads to additional "empty" result directories when running with some OpenFOAM versions, leading to inconveniences during post-processing. Please run the script `removeObsoleteSolvers.sh` to delete the additional files.

You may adjust the end time in the precice-config_*.xml, or interupt the execution earlier if you want.

This case was contributed by Derek Risseeuw (TU Delft).

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.