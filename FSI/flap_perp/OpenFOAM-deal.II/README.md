# Tutorial for an FSI simulation of an elastic flap perpendicular to a channel flow

This tutorial is described in the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-deal.II-and-OpenFOAM).

It is known to work with OpenFOAM 4.1, 5.0, v1712, and v1806, but it should also work with newer versions. Have a look into our [Notes on OpenFOAM](https://github.com/precice/openfoam-adapter/wiki/Notes-on-OpenFOAM).

You may run the coupled simulation in serial using the script `Allrun` or (OpenFOAM) in parallel with `Allrun -parallel`. The output of each step will be redirected to log files. You can cleanup the simulation using `Allclean`.

If you prefer to run the two simulations in two different terminals and watch their output on the screen, use the (simpler) scripts `runFluid` (or `runFluid -parallel`).

Before starting the deal.II program, it needs to be compiled and copied in this case directory. Information about building is collected in the [deal.II wiki](https://github.com/precice/dealii-adapter/wiki/Building). You can use the following command to run e.g. the linear deal.II solver afterwards:
```
./runSolid -linear
```

There is an [open issue](https://github.com/precice/openfoam-adapter/issues/26) that leads to additional "empty" result directories when running with some OpenFOAM versions, leading to inconveniences during post-processing. Please run the script `removeObsoleteSolvers.sh` to delete the additional files.

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
