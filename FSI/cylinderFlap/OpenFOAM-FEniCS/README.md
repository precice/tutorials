This is a tutorial for the cylinder flap case using FEniCS as a structure solver and OpenFOAM for the Fluid.
To run it you need:
- preCICE
- FEniCS
- FEniCS Adapter for FSI
- OpenFOAM + OpenFOAM Adapter (see comments on different OpenFOAM versions below)

# Run it

You can run the tutorial from this directory by typing ```python3 Solid/cyl-flap.py``` in one Terminal and ```bash ./runFluid``` in another.

## Different OpenFOAM versions

* you are using OpenFOAM 5 or older: just use the master branch of this repository together with the master branch of the OpenFOAM Adapter
* you are using OpenFOAM 6: use the branch called `OpenFOAM6` of this repository together with the branch called `OpenFOAM6` of the OpenFOAM Adapter
* you are using something else: more version of the OpenFOAM-Adapter are provided as branches of the adapter repository. However, this tutorial has not been tested with other versions of OpenFOAM.

For more information, please refer to [this issue](https://github.com/precice/tutorials/issues/40).

# Visualize Results

To visualize the results open fluid.foam in paraview. If you want to look at the flap explicitely you can open the ```.pvd``` file in paraview.

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
