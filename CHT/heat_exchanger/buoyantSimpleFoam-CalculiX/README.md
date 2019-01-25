# Tutorial for a shell-and-tube heat exchanger, using OpenFOAM and CalculiX.

This tutorial is described in the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-CHT-with-OpenFOAM-and-CalculiX).

It works with OpenFOAM 5.0 and CalculiX 2.12, but newer minor versions should work as well.

Please run the script `Download_meshes` first.
Then you may run the `Allrun` for a serial run,
or the `Allrun_parallel` for a parallel run.

You may adjust the end time in the precice-config_*.xml, or interupt the execution earlier if you want.

Based on a case prepared with simscale.com.

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
