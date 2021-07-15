# preCICE tutorials change log

All notable changes to this repository will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Read more details in the issue [#52: Releases and versioning](https://github.com/precice/openfoam-adapter/issues/52).

## [Unreleased]

## [v202104.1.1] 2021-04-23

### Added


### Changed

- Modifying the helper tool `openfoam_remove_empty_dirs` such that it now respects also results in the compressed OpenFOAM format (76f4482).
- Syncing the post-processing functionality of the elastic-tube-1d and the respective documentation. (#209)

### Removed

## [v202104.1.0] 2021-05-02

### Added

- Adding a standard run script in each case folder which can be executed as `./run.sh`.
- Adding a standard clean script in each case folder which can be executed as `./clean<what>.sh`.
- Adding an easy-to-run tutorial called [quickstart](https://precice.org/quickstart.html).
- Adding a validated Turek-Hron FSI3 case with OpenFOAM and non-linear deal.II.

### Changed

- Moving all documentation to the [redesigned preCICE website](https://precice.org/tutorials.html).
- Creating a new directory structure for easy access and also for [contributions](https://precice.org/community-contribute-to-precice.html).
- Modifying the 2D cases to use 2D mode in preCICE and also corresponding 2D functionality in the adapters.

### Removed