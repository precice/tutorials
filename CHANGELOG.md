# preCICE tutorials changelog

All notable changes to this repository will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

<!-- markdownlint-configure-file {"MD024": { "siblings_only": true } } -->

## [v202404.0] 2024-04-16

### Added

New tutorials:

- Added new multi-scale tutorial [two-scale heat conduction](https://precice.org/tutorials-two-scale-heat-conduction.html) with Nutils [#343](https://github.com/precice/tutorials/pull/343) and DuMuX [#376](https://github.com/precice/tutorials/pull/376), using the new [Micro Manager](https://precice.org/tooling-micro-manager-overview.html).
- Added new [flow around controlled moving cylinder](https://precice.org/tutorials-flow-around-controlled-moving-cylinder.html) tutorial, using the new [FMI runner](https://precice.org/tooling-fmi-runner.html) [#474](https://github.com/precice/tutorials/pull/474).
- Added new two-phase FSI tutorial [breaking dam with flexible pillar 2D](https://precice.org/tutorials-breaking-dam-2d.html) with OpenFOAM (interFoam) and CalculiX [#279](https://github.com/precice/tutorials/pull/279).
- Added new flow coupling tutorials [#326](https://github.com/precice/tutorials/pull/326):
  - [Partitioned flow over a backwards-facing step](https://precice.org/tutorials-partitioned-backwards-facing-step.html)
  - [Partitioned flow over a heated plate](https://precice.org/tutorials-flow-over-heated-plate-partitioned-flow.html)
  - [Partitioned pipe: two-phase](https://precice.org/tutorials-partitioned-pipe-two-phase.html) ([#418](https://github.com/precice/tutorials/pull/418))
- Added new volume coupling cases/tutorials with OpenFOAM:
  - Added an OpenFOAM case in the [channel transport](https://precice.org/tutorials-channel-transport.html) tutorial [#315](https://github.com/precice/tutorials/pull/315).
  - Added a new tutorial [volume-coupled flow](https://precice.org/tutorials-volume-coupled-flow.html) [#350](https://github.com/precice/tutorials/pull/350).
- Added overlapping Schwarz (Dirichlet-Dirichlet coupling) variants for the [oscillator](https://precice.org/tutorials-oscillator-overlap.html) [#391](https://github.com/precice/tutorials/pull/391) and for the [partitioned heat conduction](https://precice.org/tutorials-partitioned-heat-conduction-overlap.html) tutorials.

New solver options in existing tutorials:

- Added FMU participant options for the [oscillator](https://precice.org/tutorials-oscillator.html) tutorial [#466](https://github.com/precice/tutorials/pull/466).
- Added an SU2 case in the [flow over a heated plate](https://precice.org/tutorials-flow-over-heated-plate.html) tutorial [da7a149](https://github.com/precice/tutorials/commit/da7a1494f5c36b4ef509daf2a43bfee42fb32d9d).
- Added new solver options in the [perpendicular flap](https://precice.org/tutorials-perpendicular-flap.html) tutorial:
  - `solid-nutils` [#433](https://github.com/precice/tutorials/pull/433)
  - `fluid-fake` [#472](https://github.com/precice/tutorials/pull/472), only meant for debugging of solid participants
- Added a Rust-based solver in the [elastic tube 1D](https://precice.org/tutorials-elastic-tube-1d.html) tutorial [#435](https://github.com/precice/tutorials/pull/435).

General and documentation-related changes:

- Added automatic logging in the run scripts [#479](https://github.com/precice/tutorials/pull/479).
- Added a `.gitignore` and extended the clean-up scripts [#477](https://github.com/precice/tutorials/pull/477).
- Added visualizations of the preCICE configuration file in every tutorial [#514](https://github.com/precice/tutorials/pull/514).
- Added a list of exact Debian package links in the [quickstart](https://precice.org/quickstart.html) tutorial, to reduce cases of users installing the wrong package for their system. Now referring users to releases instead of the master branches of the Git repositories.

Developer-facing changes:

- Added on-demand system regression tests (multiple pull requests, starting from [#347](https://github.com/precice/tutorials/pull/347)).
- Added reference results, in the context of the system regression tests, hosted in an external Git LFS server [#419](https://github.com/precice/tutorials/pull/419).
- Added a pre-commit hook for preCICE configuration file formatting, and more [1d22c1](https://github.com/precice/tutorials/commit/1d22c1f61d7b13624973408c4bda7031b69adb5b), [#478](https://github.com/precice/tutorials/pull/478).

### Changed

General updates:

- Updated all example codes and configuration files to preCICE v3 (see the [porting guide](https://precice.org/couple-your-code-porting-v2-3.html)). preCICE v2 is not supported anymore. This has happened in various pull requests.
- Updated the SU2 configuration in the [perpendicular flap](https://precice.org/tutorials-perpendicular-flap.html) tutorial, for the updated SU2 adapter [da7a149](https://github.com/precice/tutorials/commit/da7a1494f5c36b4ef509daf2a43bfee42fb32d9d).
- Updated the suggested OpenFOAM version in the [Quickstart](https://precice.org/quickstart.html) to OpenFOAM v2312.
- Updated the documentation images for several tutorials, including [perpendicular flap](https://precice.org/tutorials-perpendicular-flap.html) ([#507](https://github.com/precice/tutorials/pull/507), now including reference watchpoint files), [heat exchanger: simplified](https://precice.org/tutorials-heat-exchanger-simplified.html) ([#327](https://github.com/precice/tutorials/pull/327) and [#513](https://github.com/precice/tutorials/pull/513)), [elastic tube 3d](https://precice.org/tutorials-elastic-tube-3d.html) ([#509](https://github.com/precice/tutorials/pull/509)), [multiple perpendicular flaps](https://precice.org/tutorials-multiple-perpendicular-flaps.html) ([#511](https://github.com/precice/tutorials/pull/511)), and more.

Restructured:

- Modified the following cases to fit the directory structure defined in the now extended [contributing guidelines](https://precice.org/community-contribute-to-precice.html#contributing-tutorials) ([#461](https://github.com/precice/tutorials/issues/461)):
  - [Oscillator](https://precice.org/tutorials-oscillator.html)
  - [Partitioned heat conduction](https://precice.org/tutorials-partitioned-heat-conduction.html)
  - [Partitioned heat condiction: Complex](https://precice.org/tutorials-partitioned-heat-conduction-complex.html)
  - [Partitioned heat condiction: Direct mesh access](https://precice.org/tutorials-partitioned-heat-conduction-direct.html)
  - [Volume-coupled diffusion](https://precice.org/tutorials-volume-coupled-diffusion.html)

Dependency and workflow updates:

- Modified the [Turek-Hron FSI3](https://precice.org/tutorials-turek-hron-fsi3.html) and the [partitioned heat conduction](https://precice.org/tutorials-partitioned-heat-conduction.html) OpenFOAM cases to compute the parabolic inlet profiles using a coded boundary condition, dropping the dependency on groovyBC / swak4Foam [#428](https://github.com/precice/tutorials/pull/428) (added via [a688fa](https://github.com/precice/tutorials/commit/a688fa7db044efbb72ddab7dc0bece522a5ff1e5)).
- Modified the run scripts so that Nutils is now installed automatically in a virtual environment in all related tutorials [#439](https://github.com/precice/tutorials/pull/439). Similarly for DUNE-FEM [#470](https://github.com/precice/tutorials/pull/470).
- Modified the run scripts so that C++ solvers in the elastic-tube-1d tutorial are now built automatically [#330](https://github.com/precice/tutorials/pull/330).

Improvements and bug fixes:

- Modified all [partitioned heat conduction](https://precice.org/tutorials-partitioned-heat-conduction.html) tutorials (including the basic, complex, and direct mesh access variants):
  - Modified the configuration file to use the waveform iteration (time interpolation) feature [#281](https://github.com/precice/tutorials/pull/281).
  - Modified the preCICE configuration to use a (simpler) nearest-neighbor mapping [#382](https://github.com/precice/tutorials/pull/382).
  - Modified all cases and solvers to use the same value of `beta` from 1.3 to 1.2 [#379](https://github.com/precice/tutorials/pull/379).
  - Modified the FEniCS solvers to apply boundary conditions at the right time [#383](https://github.com/precice/tutorials/pull/383).
  - Modified the FEniCS solver of the basic case to use higher-order implicit Runge-Kutta methods [#415](https://github.com/precice/tutorials/pull/415).
  - Modified the OpenFOAM solver to output error estimations [#449](https://github.com/precice/tutorials/pull/449).
- Modified the [partitioned pipe](https://precice.org/tutorials-partitioned-pipe.html) tutorial:
  - Modified the OpenFOAM cases to use the new custom boundary conditions [#326](https://github.com/precice/tutorials/pull/326).
  - Removed the unused `PressureGradient` and `VelocityGradient` coupling data [#384](https://github.com/precice/tutorials/pull/384).
- Modified the [oscillator](https://precice.org/tutorials-oscillator.html) tutorial to use higher-order time stepping schemes.
- Modified the [elastic tube 1d](https://precice.org/tutorials-elastic-tube-1d.html) tutorial to export and plot watchpoints [#438](https://github.com/precice/tutorials/pull/438).
- Modified the configuration of the [flow over a heated plate: two meshes](https://precice.org/tutorials-flow-over-heated-plate-two-meshes.html) tutorial:
  - The CalculiX domain now has the same thickness as the OpenFOAM domain [#487](https://github.com/precice/tutorials/pull/487).
  - Removed an unused mesh [#390](https://github.com/precice/tutorials/pull/390).

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
