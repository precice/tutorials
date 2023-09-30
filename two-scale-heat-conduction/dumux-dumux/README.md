# Two-scale heat conduction with DuMu<sup>x</sup> and preCICE

This repository showcases two-scale coupling of [DuMu<sup>x</sup>](https://dumux.org/) simulations with [preCICE](https://www.precice.org/), using a two-scale coupled heat conduction problem in 2D as an example. The DuMu<sup>x</sup> simulations can also be coupled with [their respective counterparts in Nutils](https://github.com/IshaanDesai/coupled-heat-conduction).

## Setup

Being an example, this repository is sensitive to versioning, so make sure to get the correct versions of the involved software packages. **TLDR**s might be sufficient for quick setup, but they ignore all dependencies and will almost certainly fail or create conflicts for any system without previous configuration. Consulting the full setup instructions for each software component is very much recommended.

After following the steps below, your repository structure should look like this:

```bash
|-- dumux
    |-- dumux
        |-- dumux-adapter
            |--examples
               |--macro-micro
               ...
        |-- dumux-phasefield
        |-- dune-common
        |-- dune-geometry
        |-- dune-grid
        |-- dune-istl
        |-- dune-localfunctions
        |-- dune-spgrid
```

## Dependencies

Apart from DuMu<sup>x</sup> and preCICE, this project relies on three other software components: the [DuMu<sup>x</sup>-preCICE adapter](https://github.com/precice/dumux-adapter/), the [Micro Manager](https://github.com/precice/micro-manager) and the (currently still private) DuMuX module [dumux-phasefield](https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield/).

### Micro Manager

**Version**: `develop` branch

Follow the [installation instructions for the Micro Manager](https://github.com/precice/micro-manager/blob/main/README.md). This also lists the necessary dependencies such as [preCICE](https://www.precice.org/).

**Note**: Currently (at the time of precice-micro-manager v.0.2.1), this code relies on the [develop](https://github.com/precice/micro-manager/tree/develop) branch, so make sure to follow the steps for [manual installation](https://github.com/precice/micro-manager/tree/main#option-2-clone-this-repository-and-install-manually) and use the `develop` branch.

**TLDR**:

```bash
git clone -b develop https://github.com/precice/micro-manager
cd micro-manager
pip install --user .
```

### DuMu<sup>x</sup>

**Version**: DuMu<sup>x</sup> 3.6, DUNE 2.8

Follow the [installation instructions for DuMuX](https://dumux.org/installation/). Using the install script is recommended, but check that it contains the correct versions. This creates the following directory structure:

```bash
|-- dumux
    |-- dumux
    |-- dune-common
    |-- dune-geometry
    |-- dune-grid
    |-- dune-istl
    |-- dune-localfunctions

```

In the following, `dumux` always refers to the *outer* dumux folder.

**TLDR**:

Download the [install script](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/bin/installdumux.py) and run it with `python installdumux.py`.

### DuMu<sup>x</sup>-preCICE adapter

**Version**: precice/dumux-adapter v.1.0.0

Follow the [installation instructions for the DuMu<sup>x</sup> adapter](https://github.com/precice/dumux-adapter/).

**TLDR**:

```bash
cd dumux
git clone -b v1.0.0 https://github.com/precice/dumux-adapter.git
./dune-common/bin/dunecontrol --only=dumux-adapter/dumux-precice all
```

### Additional DuMu<sup>x</sup> Modules

**Version**: dune-spgrid 2.6, dumux-phasefield `cell_problems` branch

Steps to download and build the modules:

1. Clone `dumux-phasefield` and `dune-spgrid` into the dumux folder.

```bash
cd dumux
git clone -b cell_problems https://git.iws.uni-stuttgart.de/dumux-appl/dumux-phasefield/
python3 dumux/bin/installexternal.py spgrid
```

**Note**: `dumux-phasefield` is still private!

2. Clean the CMake cache.

```bash
./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt
```

3. Reconfigure and build DuMu<sup>x</sup> via dunecontrol.

```bash
./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all
```

## Setup of the macro-micro example

**Version**: macro-micro `main` branch , pybind11 v2.10

1. Navigate into `dumux-adapter/examples/` and clone this repo.

```bash
cd dumux/dumux-adapter/examples/
git clone https://github.com/HelenaKschidock/macro-micro
```

2. Add the line `add_subdirectory(macro-micro)` to `dumux-adapter/examples/CMakeLists.txt`. This integrates the example into the overall cmake build  structure.

3. Add `dumux-phasefield` to `dumux-adapter/dune.module` as a required build dependency. The line should then look like this:

```
# Required build dependencies
Depends: dumux (>=3.2) dumux-phasefield
```

4. Navigate into the `micro-heat` directory and clone and build pybind11 there. For more information, see also [https://github.com/pybind/cmake_example](https://github.com/pybind/cmake_examples).

```bash
cd macro-micro/micro-heat
git clone -b v2.10 https://github.com/pybind/pybind11
cd pybind11 
pip install --user .
```

5. Rebuild the DuMu<sup>x</sup>-preCICE adapter.

```bash
cd dumux
./dune-common/bin/dunecontrol --only=dumux-precice all
```

### Building the macro and micro executables

To build the macro and micro executables, simply rebuild the dumux-adapter using dunecontrol:

```bash
cd dumux 
./dune-common/bin/dunecontrol --only=dumux-precice all
```

The executables are built into `dumux/dumux-adapter/build-cmake/examples/macro-micro`.

In case you want to rebuild only the macro executable, this can also be done directly via:

```bash
cd dumux/dumux-adapter/build-cmake/examples/macro-micro/macro-heat
make test_macro_heat
```

## Running The Coupled Simulation

0. Open two terminals.

1. In one of them, run the macro simulation via

```bash
cd dumux/dumux-adapter/build-cmake/examples/macro-micro/macro-heat 
./test_macro_heat
```

2. In the other shell, run the micro simulation via

```bash
cd dumux/dumux-adapter/build-cmake/examples/macro-micro/micro-heat
python3 run_micro_manager.py 
```

To run the micro code in parallel, use `mpirun -n <NThreads> python3 run-micro-manager.py` instead. **(not tested yet)**

## Coupling with Nutils

For the same problem, similar macro and micro simulations have also been [implemented in Nutils](https://github.com/IshaanDesai/coupled-heat-conduction), which can be coupled with their DuMu<sup>x</sup> counterparts. To do so:

### Setup the Nutils macro and micro case files

1. [Install Nutils](https://nutils.org/install.html).

2. Clone the coupled-heat-conduction git into your macro-micro directory.

```bash
cd dumux-adapter/examples/macro-micro
git clone https://github.com/IshaanDesai/coupled-heat-conduction
```

3. Add the line `add_subdirectory(coupled-heat-conduction)` to `macro-micro/CMakeLists.txt` to link it to CMake.

4. To make the coupled-heat-conduction code available in your build directory, create a new empty `CMakeLists.txt` file within your source `coupled-heat-conduction` directory and add the following code:

```cmake
add_custom_target(copy_nutils_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
```

This simply copies it to the build directory and time you execute the `copy_nutils_all` command or anytime you rebuild the entire `dumux-precice` module.

5. Modify the `coupled-heat-conduction/precice-config.xml` file such that it reads

```xml
 <m2n:sockets from="Micro-Manager" to="Macro-heat" network="lo" exchange-directory="../"/>
```

This allows the macro and micro simulations to find each other.

6. Rebuild the adapter.

```bash
dunecontrol --only=dumux-precice all
```

### Coupling the Nutils simulations

You can couple the Nutils simulations with each other or with their counterpart in DuMu<sup>x</sup>, as long as one participant is the macro simulation and the other the micro simulation/micro manager. To do so, simply execute your chosen participants in separate shells.

* Run the Nutils macro simulation via

```bash
cd build-cmake/examples/macro-micro/coupled-heat-conduction
python3 macro-heat.py
```

* Run the micro simulation via

```bash
cd build-cmake/examples/macro-micro/coupled-heat-conduction
python3 run-micro-manager.py
```

or in parallel via

```bash
mpirun -n <num_procs> python3 run-micro-manager.py
```

## TBC

## Note about the Macro parameters

In order to modify the simulation parameters, modify `params.input`.
The parameters are currently tuned to the dimensionless micro-simulation. When modifying these, note the following:

* The domain is defined by its `LowerLeft` and `UpperRight` corners and split into `[nx ny] Cells`.
* Set `RunWithCoupling` to `false` to run the macro simulation independently or to `true` to run the coupled simulation. The independent simulation uses `SolidThermalConductivity` and `DefaultPorosity`, respectively, as its conductivity tensor and porosity.
* Set `BcTypeLeft/Right/Top/Bottom` to `dirichlet` or `neumann`, and the corresponding value in `BcLeft/Right/Top/Bottom`. Note that Neumann boundary conditions are currently hardcoded to 0.0.
* To run the hardcoded asymmetric case with a heat source in the bottom left, set `UseHeatSourceBottomLeft` to `true`, else to `false`. The corresponding value can be set via `HeatSourceBottomLeft`.
