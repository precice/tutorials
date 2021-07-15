# preCICE tutorials change log

All notable changes to this repository will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Read more details in the issue [#52: Releases and versioning](https://github.com/precice/openfoam-adapter/issues/52).

<!-- markdownlint-configure-file {"MD024": { "siblings_only": true } } -->

## [Unreleased]

## [v202104.1.1] 2021-04-23

### Changed

- Modified the helper tool `openfoam_remove_empty_dirs` such that it also respects results in the compressed OpenFOAM format (76f4482).
- Synced the post-processing functionality of the elastic-tube-1d and the respective documentation. (#209)

## [v202104.1.0] 2021-05-02

### Added

- Added a change log for this project as a file named `CHANGELOG.md`.
- Added a standard run script in each case folder which can be executed as `./run.sh`.
- Added a standard clean script in each case folder which can be executed as `./clean<what>.sh`.
- Added an easy-to-run tutorial called [quickstart](https://precice.org/quickstart.html).
- Added a validated Turek-Hron FSI3 case with OpenFOAM and non-linear deal.II.

### Changed

- Moved all documentation to the [redesigned preCICE website](https://precice.org/tutorials.html).
- Created a new directory structure for easy access and also for [contributions](https://precice.org/community-contribute-to-precice.html).
- Modified the 2D cases to use 2D mode in preCICE and also corresponding 2D functionality in the adapters.
