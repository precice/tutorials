# preCICE tutorials change log

All notable changes to this repository will be documented in this file.

## 202104.1.1

- The modified helper tool `openfoam_remove_empty_dirs` now respects also results in the compressed OpenFOAM format (76f4482).
- Syncing the post-processing functionality of the elastic-tube-1d and the respective documentation. (#209)

## 202104.1.0

- **All documentation on the website:** Find everything directly in our [redesigned preCICE website](https://precice.org/tutorials.html) (rendered from the `README.md` files here, so you can also check the basic information without internet connection).
- **New directory structure:** Read more in the [contributing guidelines](https://precice.org/community-contribute-to-precice.html).
- **Self-contained cases:** The files for each case are inside its own directory, e.g. `fluid-openfoam`. The exchange directory is set to `..`, next to `precice-config.xml`.
- **Arbitrary combinations of solvers:** Run OpenFOAM-CalculiX, OpenFOAM-deal.II, SU2-deal.II
- **A standard way to run each case:** Go to the case folder and execute `./run.sh`.
- **A standard way to clean:** A `clean<what>.sh` script in each directory cleans up everything from this level and deeper.
- **An easy-to-run first tutorial:** We created a [quickstart](https://precice.org/quickstart.html) tutorial.
- **Better 2D cases:** Two-dimensional cases such as the (now rotated) perpendicular flap are using the 2D interface of preCICE and every involved adapter supports a 2D mode.
- **A validated Turek-Hron FSI3:** Validated with OpenFOAM and non-linear deal.II.

