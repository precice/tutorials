# Tutorial for an CHT simulation of a flow over a heated plate using OpenFOAM and Nutils

## Case setup: flow over heated plate

The setup for this tutorial is identical to the [flow over heated plate using OpenFOAM](https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate#case-setup) with the only exception that the width of the scenario in z-direction is `0.05` instead of `0.4`. 

## Dependencies

### Nutils 

[Nutils](http://www.nutils.org/en/latest/) is an open-source Python programming finite element library, developed by [Evalf Computing](http://evalf.com/). 

Clone and install via pip:

```
$ git clone https://github.com/nutils/nutils.git
$ python3 -m pip install --user --editable nutils
```

For faster computations, you can optionally install `mkl`:

```
$ pip3 install mkl
```

### Other dependencies

For running this tutorial, you further have to install

* **preCICE**, see [preCICE wiki](https://github.com/precice/precice/wiki/Building).
* **Python bindings**, see [bindings readme](https://github.com/precice/precice/tree/develop/src/precice/bindings/python)
* **OpenFOAM**, see [Notes on OpenFOAM](https://github.com/precice/openfoam-adapter/wiki/Notes-on-OpenFOAM).
* **OpenFOAM adapter**, see [OpenFOAM adapter wiki](https://github.com/precice/openfoam-adapter/wiki/Building). If you have problems compiling, see the [troubleshooting section](https://github.com/precice/precice/wiki/CHT-with-OpenFOAM-and-FEniCS#troubleshooting) below.

### Testing your installation

* **OpenFOAM and OpenFOAM adapter:** To make sure that everything is working properly, you should run the following tutorial case: https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate.
* **Nutils:** To make sure that Nutils is working properly, you should run at least one of the [Nutils examples](http://www.nutils.org/en/latest/examples/).

## Run the tutorial

Open two terminals at the root of this tutorial.

Terminal 1:
```
$ cd OpenFOAM
$ ./runFOAM
```

Terminal 2:
```
$ cd Nutils
$ python3 cht.py
```

Alternatively, you can also directly use the `Allrun` script in one terminal. 

### Visualization

Both solvers, OpenFOAM and Nutils, create vtk output that you can, for example, load in Paraview. 

After 100 timesteps with `dt=0.01`: 

![Visualization of the temperature](https://raw.githubusercontent.com/wiki/precice/precice/images/CHT_OpenFOAM_Nutils.png)

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
