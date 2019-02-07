# Experimental Setup

We want to compare the flow over heated plate scenario for two different software setups. We also consider the analytic solution for comparison. We only consider the physical setup with **k=1**, **lambda=0.25**, **Re=500**, **Pr=0.01**, **T_c = 310**, **T_inf=310**. For explanation of the parameters and reference solution, refer to [1]. This specific parameter set has been evaluated in [1] Fig.8 a), [2], [3] and [4].

## 1) OF-OF

We use OpenFOAM for simulation of the fluid and OpenFOAM for simulation of the solid

* source code and running, see [here](LINK!!!)
* preCICE config being used [`precice-config_serial_FE_OF.xml`](LINK!!!)
* data stored in `out_OF-OF.csv`

## 2) FE-OF

Using OpenFOAM for simulation of the fluid and FEniCS for simulation of the solid.

* source code and running, see [here](LINK!!!)
* preCICE config being used [`precice-config_serial_FE_OF.xml`](LINK!!!)
* data stored in `out_FE-OF.csv`

## 3) analytic solution

For derivation see [1].

## Comments

The data `*.out` is obtained from the vtk output of the solid part of the simulation. We use paraview to exporting the data of *plot over line* along the coupling boundary.

The simulation has to be run for a larger amount of time (approx. 50s simulation time, see `precice-config`) in order to reach steady state. **For the sake of speeding up the computation, the timestep size has been increased for the FEniCS solver!**

## Postprocessing

Run `python3 plotParaviewOut.py` to obtain the plot with the results.
 
## Results

See `comparison.png`

## References

[1] Vynnycky Paper PUT PROPER REFERENCE!
[2] Cheung Thesis PUT PROPER REFERENCE!
[3] Reiser Thesis PUT PROPER REFERENCE!
[4] Chourdak Thesis PUT PROPER REFERENCE!
