# Two Perpendicular Flap

In the following tutorial we model a fluid flowing through a channel. 
Two solid, elastic flaps are fixed to the floor of this channel.
The flaps oscillate due to the fluid pressure building up on its surface. This case is an example for multi-coupling: a fluid and two solids are coupled together using a fully-implicit multi-coupling scheme.

## Case Setup

The setup is NOT YET shown schematically here:

[[images/setupDrawing.png]]

The simulated flow domain is 6 units long (x) and 4 units tall (z). The flaps are clamped at the bottom (z=0) and they are 1 unit tall (z), 0.1 units long (x), and 0.3 units wide (y). Being located at x=-1 and x=1, the flaps split the domain into three equal parts. 

This is a quasi-2D simulation, with one cell in width and empty conditions are imposed on the front and back. 

The inflow velocity is 10 m/s (uniform) on the left boundary.
At the outlet, pressure is set to zero and velocity to `zeroGradient`.
The top, bottom and flap are walls with the `noslip` condition. 


## Running the Simulation
SECTION TO BE UPDATE
1. Preparation:
   To run the coupled simulation, copy the deal.II executable `linear_elasticity` or `nonlinear_elasticity` into the `Solid` folder.           For OpenFOAM: The name of your solver might differ, depending on your OpenFOAM version. Have a look in the `Fluid/system/controlDict` file and set the appropriate solver name.
2. Starting:
   We prepared a shell script to start the `Fluid` participant. Run
   ```
   ./runFluid
   ```
   to start OpenFOAM in serial or
   ```
   ./runFluid -parallel
   ```
   for a parallel run.
   The Solid participant is only designed for serial runs. Therefore, open another terminal, navigate to the simulation directory and execute the deal.II binary file e.g. by:
   ```
   ./linear_elasticity ./Solid/linear_elasticity.prm
   ```
   or simply `./runSolid -linear`. Similarly for the non-linear case.
   It is important to start both solvers in the same directory.
## Postprocessing
SECTION TO BE UPDATE

After the simulation has finished, you can visualize your results using e.g. ParaView. Fluid results are in the OpenFOAM format and you may load the `Fluid.foam` file. Solid results are in VTK format and located in the `dealii_output` directory. If you want to plot both domains with ParaView, keep in mind that the deal.II solver we presented here writes results every few timesteps, while the OpenFOAM solver writes in reference to simulated time. For this reason, make sure that you use compatible write intervals. You may also need to convert the OpenFOAM results to VTK (with the command `foamToVTK`).

Apart from the usual result files, a watchpoint has been defined in the `precice-config.xml` file, which allows one to track coupling data at a specific point at the coupling interface. We chose the y-displacement of the flap tip:

![flap_tip](https://user-images.githubusercontent.com/33414590/58786669-62eb3780-85e8-11e9-81a0-5432c3e251c3.png)

Besides the obtained results of our coupled deal.II-OpenFOAM case, the graph shows reference data, which have been calculated with a geometrical nonlinear structural model. Here, it is easy to see, that our linear model is not sufficient for large displacements. Therefore, it is recommended to use the nonlinear solver.

# References

[1]  S. Turek, J. Hron, M. Madlik, M. Razzaq, H. Wobker, and J. Acker. Numerical simulation and benchmarking of a monolithic multigrid solver for fluid-structure interaction problems with application to hemodynamics. In H.-J. Bungartz, M. Mehl, and M. Schäfer, editors, Fluid Structure Interaction II: Modelling, Simulation, Optimization, page 432. Springer Berlin Heidelberg, 2010.
