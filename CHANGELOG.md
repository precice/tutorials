# preCICE tutorials changelog

All notable changes to this repository will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

<!-- markdownlint-configure-file {"MD024": { "siblings_only": true } } -->

## [v202211.0] 2022-11-21

### Added

- Added an [ASTE tutorial](https://precice.org/tutorials-aste-turbine.html) [#244](https://github.com/precice/tutorials/pull/244).
- Added an [oscillator tutorial](https://precice.org/tutorials-oscillator.html) with Python [#297](https://github.com/precice/tutorials/pull/297).
- Added a uni-directional volume-coupled [channel transport tutorial](https://precice.org/tutorials-channel-transport.html) [#269](https://github.com/precice/tutorials/pull/269).
- Added a volume-coupled [channel transport with chemical reactions tutorial](https://precice.org/tutorials-channel-transport-reaction.html) [#278](https://github.com/precice/tutorials/pull/278).
- Added a [partitioned heat conduction using direct mesh access tutorial](https://precice.org/tutorials-partitioned-heat-conduction-direct.html) [#299](https://github.com/precice/tutorials/pull/299).
- Added a [simplified heat exchanger tutorial](https://precice.org/tutorials-heat-exchanger-simplified.html) [#301](https://github.com/precice/tutorials/pull/301).
- Added a solid-openfoam and a solid-solids4foam case for the perpendicular-flap tutorial [#286](https://github.com/precice/tutorials/pull/286).
- Added a dune-fem case for the flow-over-heated-plate tutorial [#274](https://github.com/precice/tutorials/pull/274).
- Added a CalculiX case for the flow-over-heated-plate tutorial [#271](https://github.com/precice/tutorials/pull/271).
- Added a modal dynamic simulation option for the CalculiX case of the perpendicular-flap tutorial [#284](https://github.com/precice/tutorials/pull/284).
- Added a CI workflow to automatically update the website every time there are new documentation changes [#267](https://github.com/precice/tutorials/pull/267)
- Added a gitignore for code and binaries [#290](https://github.com/precice/tutorials/pull/290).
- Added more documentation [#264](https://github.com/precice/tutorials/pull/264), [#265](https://github.com/precice/tutorials/pull/265), [#266](https://github.com/precice/tutorials/pull/266) and fixed some typos [#285](https://github.com/precice/tutorials/pull/285).

### Changed

- Adapted the multiple-perpendicular-flaps tutorial to a [new naming convention](https://precice.org/community-contribute-to-precice.html#contributing-tutorials) for tutorials with multiple participants [#303](https://github.com/precice/tutorials/pull/303).
- Adapted the material properties of the FEniCS case in elastic-tube-3d to match the CalculiX case results [#259](https://github.com/precice/tutorials/issues/259).
- Adapted the solid case of the Quickstart tutorial to work with the Intel oneAPI C++ compiler [#291](https://github.com/precice/tutorials/pull/291).
- Made the Nutils case of the partitioned heat conduction tutorial compatible with Nutils versions 6 to 8 [#295](https://github.com/precice/tutorials/pull/295), [#298](https://github.com/precice/tutorials/pull/298).

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
- Removed unnecessary (wrong) read statement in `elastic-tube-1d` [#232](https://github.com/precice/tutorials/pull/232).
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
