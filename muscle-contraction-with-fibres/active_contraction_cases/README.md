# Active contraction cases

This folder contains the final muscle contraction cases.

## Cases
The following cases are provided:
- [muscle_only](muscle_only): A muscle without any fibers solved using OpenDiHus MuscleContractionSolver.
- [fibers_only](fibers_only): Activation of multiple fibers solved using OpenDiHus FastMonodomainSolver.
- [muscle_contraction](muscle_contraction): Full muscle simulation with OpenDiHu using OpenDiHus internal coupling.
- [precice_contraction](precice_contraction): Full muscle simulation with OpenDiHu using external preCICE coupling.
- [febio_contraction](febio_contraction): Couples FEBios mechanics solver with OpenDiHus FastMonodomain solver using the preCICE adapter.

## Running the OpenDiHu cases
To run an OpenDiHu case, you first have to compile it using
```bash
opendihu/scripts/shortcuts/sr.sh
```
After this you can navigate to *build_release* and run
```bash
./program_name ../settings_program_name.py config_file
```
where *config_file* is one of the following:

### variables.py
Default configuration file. Uses a 3x3x12 mesh and 0.1ms timesteps for the mechanics part.
The fibers have a resolution of 10x10 and use 100 points per fiber.
The mesh is fixed to the z=0 plane using Dirichlet boundary conditions.
The simulation runs for 30ms. 

### neumann.py
The same as *variables.py*, but includes Neumann pulling force in z direction.
Supposed to be used with the *muscle_only* case.

### mapping.py
Similar to *variables.py*, but with a 4x4 fiber resolution.
With this case the internal OpenDiHu mapping method does not work properly.
This means that the results of *precice_contraction* and *muscle_contraction* will differ significantly.

### case_zN_Tms.py
Similar to *variables.py*, but uses N elements in z direction and a timestep width of T ms instead.
This can be used to compare the results at different resolutions.

## Running the FEBio cases
The FEBio cases are located in *febio_contraction/muscle*.
You can change the *run.sh* to point to your FEBio4 executable.
The case *muscle.feb* is supposed to be used with the default *variables.py* file on the OpenDiHu side,
and the cases *muscle_zN_Tms.feb* are supposed to be used with the corresponding *case_zN_Tms.py*:
```bash
BFP_CONFIG="../../variables/precice-config-Tms.xml" ./run.sh muscle_zN_Tms.feb
./fibers ../settings_fibers.py case_zN_Tms.py
```
Note that you need to set the *BFP_CONFIG* environment variable if you want the preCICE adapter to use the correct configuration.

