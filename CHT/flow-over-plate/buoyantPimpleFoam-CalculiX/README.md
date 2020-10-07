# Tutorial for an CHT simulation of a flow over a heated plate using OpenFOAM and CalculiX

## Case setup: flow over heated plate

The setup for this tutorial is identical to the [flow over heated plate using OpenFOAM](https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate#case-setup) with the only exception that the width of the scenario in z-direction is `0.05` instead of `0.4`. 



### Other dependencies

For running this tutorial, you further have to install

* **preCICE**, see [preCICE wiki](https://github.com/precice/precice/wiki/Building).
* **Python bindings**, see [`precice/python-bindings`](https://github.com/precice/python-bindings)
* **OpenFOAM**, see [Notes on OpenFOAM](https://github.com/precice/openfoam-adapter/wiki/Notes-on-OpenFOAM).
* **OpenFOAM adapter**, see [OpenFOAM adapter wiki](https://github.com/precice/openfoam-adapter/wiki/Building). If you have problems compiling, see the [troubleshooting section](https://github.com/precice/precice/wiki/CHT-with-OpenFOAM-and-FEniCS#troubleshooting) below.
* **CalculiX adapter**, see [CalculiX adapter wiki](https://github.com/precice/calculix-adapter/wiki)


### Testing your installation

* **OpenFOAM and OpenFOAM adapter:** To make sure that everything is working properly, you should run the [similar OpenFOAM-OpenFOAM tutorial case](https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate).

### Generating CacluliX files
After installing the CalculiX adapter, CalculiX on it's own should also be available in the system. If not, an easy to use precomplied [tool](http://www.calculixforwin.com/) may be used to produce the required files (the commands remain same, use the pre-processing option).
The files for CalculiX are already generated from the `Mesh_Coarse_OUT.inp` by using the command in the Solid directory.

Terminal  :
```
Solid/$ cgx -c Mesh_Coarse_OUT.inp
```

This will launch a cgx window with geometry on it. Stay on this window and execute the following commands:

* To show all the groups available:

```
prnt se
```

* To generate the mesh and node files:

```
send all abq
send fix1 abq nam
send interface abq nam
```

* To generate the surface file

```
send sinterface abq sur
```

* To generate the files for flux and convection film

```
send sinterface abq dflux 0.0
send sinterface abq film 300 4000
```

## Run the tutorial

Open two terminals at the root of this tutorial.

Terminal 1:
```
$ ./runFluid
```

Terminal 2:
```
$ ./runSolid
```

Alternatively, you can also directly use the `Allrun` script in one terminal. 

### Visualization

OpenFOAM creates vtk output that you can, for example, load in Paraview. For CalculiX, one can use [ccx2paraview](https://github.com/calculix/ccx2paraview) to convert `.frd` to `.vtk` or `.vtu`.


## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
