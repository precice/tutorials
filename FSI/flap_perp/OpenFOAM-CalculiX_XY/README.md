## flap_perp_Salome_XY


- In this case, all the geometries and meshes are created with Salome
- There are a few changes made to `runFluid`, `runSolid`, `Allrun` and `Allclean`. These particularly have some commands commented out so as to not delete important files.
- This case has same orientation as that of preCICE tutorial case.
- This case is oriented in the XY direction.

### To run this case

- Either run `./Allrun` in a terminal or
- run the Fluid and Solid participants in different terminals using the commands `./runFluid` and `./runSolid` respectively


### Fluid

- The geometry Salome file (`.hdf`) and the mesh file (`.unv`) are available. If by chance an error comes where `polymesh` is missing, in the Fluid directory run `ideasUnvToFoam Mesh_name.unv` from the case directory, this should create the required files
- If there is an error for frontAndBack (`empty`), modify `Fluid/constant/polymesh/boundary` file as per precice case.
- The Fluid folder is used from preCICE flap tutorial case.
- `blockMeshdict` is renamed.


### CalculiX

- The `flap.hdf` is the Salome file. `flap.unv` is the mesh file for CalculiX.
- use [`CalculiX4Caelinux` CalculiX launcher](http://calculixforwin.blogspot.com/2015/05/calculix-launcher.html) or [`unv2ccx`](https://github.com/calculix/unv2ccx/releases)
- Once the `.unv` file is converted. Copied the `flap.inp` file from preCICE tutorial case. In the `CalculiX4Caelinux`, open the `Mesh_name_OUT.inp` in pre-processor. 
- While on the `cgx` window, type `prnt se`. This will list out all the group names.
- Use following commands
```
send all abq
send fix abq nam
send interface abq nam

sys mv interface.nam interface_beam.nam
sys mv fix.nam fix1_beam.nam
```

- Modified the names of the datasets as need. Adapt `config.yml`.
- Modified the names of boundary, load etc in `flap.inp`




