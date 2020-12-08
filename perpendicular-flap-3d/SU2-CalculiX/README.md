# Tutorial for an FSI simulation of an elastic flap perpendicular to a channel flow using SU2 and CalculiX

This tutorial is described in the [preCICE wiki](https://github.com/precice/precice/wiki/FSI-tutorial).

You may run the coupled simulation in serial using the script `Allrun` or in parallel with `Allrun -parallel`. The output of each step will be redirected to log files (mainly `Fluid.log` and `Solid.log`). You can cleanup the simulation using `Allclean`.

If you prefer to run the two simulations in two different terminals and watch their output on the screen, use the (simpler) scripts `runFluid` (or `runFluid -parallel`) and `runSolid`. They also write the output to the same log files.
