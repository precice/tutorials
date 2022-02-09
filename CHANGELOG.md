# preCICE tutorials changelog

All notable changes to this repository will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

<!-- markdownlint-configure-file {"MD024": { "siblings_only": true } } -->

## [Unreleased]

## [v202202.0] 2022-02-09

### Added

- Added new volume-coupled-diffusion tutorial with FEniCS [#219](https://github.com/precice/tutorials/pull/219).
- Added OpenFOAM case to partiitoned-heat [#223](https://github.com/precice/tutorials/pull/223).
- Added DUNE case to perpendicular-flap [#239](https://github.com/precice/tutorials/pull/239).
- Added FEniCS case to elastic-tube-3d [#223](https://github.com/precice/tutorials/pull/223).
- Added this changelog, describing also changes in previous releases [#225](https://github.com/precice/tutorials/pull/225).

### Changed

- Changed C3D8 elements to C3D8I elements in perpendicular-flap solid-calculix to improve the results [#250](https://github.com/precice/tutorials/pull/250).
- Ported the `visualize.py` script of the partitioned-elastic-beam to Python 3 [#247](https://github.com/precice/tutorials/pull/247).
- Reduced the writing frequency of the partitioned-pipe OpenFOAM cases [#257](https://github.com/precice/tutorials/pull/257).
- Renamed the output directories of all FEniCS cases for consistency [#256](https://github.com/precice/tutorials/pull/257).
- Removed unnecessary (wrong) read statment in `elastic-tube-1d` [#232](https://github.com/precice/tutorials/pull/232).
- Removed unnecessary (relic) OpenFOAM parameter `nMoles` from flow-over-heated-plate cases with OpenFOAM [#234](https://github.com/precice/tutorials/pull/234).
- Removed unnecessary (relic) OpenFOAM parameter `RAS` from the `turbulenceProperties` files of all OpenFOAM cases (we model a laminar flow everywhere and this was confusing) [#258](https://github.com/precice/tutorials/pull/258).
- Removed unnecessary (relic) OpenFOAM files `RASProperties` and `couplingProperties` from elastic-tube-3d [#258](https://github.com/precice/tutorials/pull/258).
- Removed unnecessary (relic) OpenFOAM file `radiationProperties` from heat-exchanger [#258](https://github.com/precice/tutorials/pull/258).
- Removed duplicate default settings from the `fvSchemes` OpenFOAM file of the heat-exchanger [#258](https://github.com/precice/tutorials/pull/258).
- Removed unnecessary (relic) fields and inaccurate copyright notices from the headers of several OpenFOAM files [#258](https://github.com/precice/tutorials/pull/258).
- Adjusted the formatting in several OpenFOAM files [#258](https://github.com/precice/tutorials/pull/258).
- Changed the versioning scheme from `<yearmonth>.<minor>.<bugfix>` to `<yearmonth>.<bugfix>`.

## [v202104.1.1] 2021-05-02

### Changed

- Modified the helper tool `openfoam_remove_empty_dirs` such that it also respects results in the compressed OpenFOAM format (76f4482).
- Synced the post-processing functionality of the elastic-tube-1d and the respective documentation. (#209)

## [v202104.1.0] 2021-04-23

### Added

- Created a first tagged version of this repository along with a [release](https://github.com/precice/tutorials/releases/tag/v202104.1.0).
- Added a standard run script in each case folder which can be executed as `./run.sh`.
- Added a standard clean script in each case folder which can be executed as `./clean<what>.sh`.
- Added an easy-to-run tutorial called [quickstart](https://precice.org/quickstart.html).
- Added a validated Turek-Hron FSI3 case with OpenFOAM and non-linear deal.II.

### Changed

- Moved all documentation to the [redesigned preCICE website](https://precice.org/tutorials.html).
- Created a new directory structure for easy access and also for [contributions](https://precice.org/community-contribute-to-precice.html).
- Modified the 2D cases to use 2D mode in preCICE and also corresponding 2D functionality in the adapters.
