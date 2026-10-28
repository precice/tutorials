# Pulsating Flow over Flexible Flap: 3D FSI Coupling Test Case

This is a variant of the "Pulsating flow over flexible flap" test case discussed in the [OpenFOAM adapter validation paper](https://journal.openfoam.com/index.php/ofj/article/view/88/78) where the `Solid`
participant is modelled in OpenRadioss.

Additionally, the OpenFOAM adapter validation paper originally aimed to validate an OpenFOAM - deal.II multiphysics coupling scenario by recreating the `pulsating flow over a bending flexible flap` monolithic simulation presented in the [Monolithic cut finite element–based approaches for fluid-structure interaction](https://onlinelibrary.wiley.com/doi/10.1002/nme.6072) paper; in this repository the `Solid` participant configured for the coupling has been replaced by an OpenRadioss comparable one.

Given the current restrictions presented when interfacing OpenRadioss (an explicit solver) to the preCICE coupling library, the `coupling-scheme` choice implemented in this repository
fell on a `parallel-explicit` configuration over the `parallel-implicit` scheme implemented in the paper mentioned above.

![preCICE configuration visualization](images/pulsating-flow-over-flexible-flap-3d-precice-config.png)

## Mesh Configuration

The test case implemented within this repository implements what the [OpenFOAM adapter validation paper](https://journal.openfoam.com/index.php/ofj/article/view/88/78) refer to as `coarse` mesh.

This means the `Solid` participant modelled in OpenRadioss is discretized with 48 `HEX20` solid elements, arranged with a 2 x 4 x 6 cell resolution along the three spatial directions.

## How to run

1. Install OpenRadioss and set the `OPENRADIOSS_PATH` in your environment (includes the `starter/` and `engine/` directories)
2. Install the `anim_to_vtk` in `OPENRADIOSS_PATH/exec`
3. Install [OpenFOAM](https://www.openfoam.com/) the [OpenFOAM adapter](https://precice.org/adapter-openfoam-overview.html) (v1.4.0 or later)
4. Open two terminals and execute the `run.sh` scripts:

   ```sh
   # Terminal 1
   cd fluid-openfoam
   ./run.sh

   # Terminal 2
   cd solid-openradioss
   ./run.sh
   ```

   The scripts rely on tools found in `../tools/`.

## State and results

The simulation completes, and the results are as follows:

### Simulation setup at `t= 2.5 s`

<img width="1920" height="1080" alt="1_Set-up-2and5secs-1" src="https://github.com/user-attachments/assets/7e4a5dc1-ab3d-4dd6-8717-9cb1baa18a57" />
<img width="1920" height="1080" alt="2_Set-up-2and5secs-2" src="https://github.com/user-attachments/assets/44b8482e-d1d2-4d97-82dc-c1a5eadfe597" />
<img width="1920" height="1080" alt="3_Set-up-2and5secs-3" src="https://github.com/user-attachments/assets/aeb9b4ea-9db9-44e8-8a1a-c103f90cfd9f" />

### Cross-section of the computational domain at `z= 0 m` and `t= 2.5s`

![4_cross-section-2and5secs](https://github.com/user-attachments/assets/2cb79f3c-ff19-4652-a4e6-5dc313e4eb04)

### Cross-section of the computational domain at `z= 0 m`: Temporal Evolution

https://github.com/user-attachments/assets/93a2f6d9-3868-43ab-ae0d-e25be37def33

### Watchpoints post-processing

Central WatchPoint

![watchpoint_Center](https://github.com/user-attachments/assets/221b0fb6-95aa-4359-922f-29a32a7b3b56)

Left WatchPoint

![watchpoint_Left](https://github.com/user-attachments/assets/af9a918b-7451-4f9a-b67f-6a5b09005258)
