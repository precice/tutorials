# Reference watchpoints

Produced by Gerasimos Chourdakis on a local system with a setup close to the
[preCICE Distribution](https://precice.org/installation-distribution.html) v2404,
as testing overlaps with the release process:

- Ubuntu 22.04
- preCICE v3.1.0 release candidate from March 28 (commit `3a1e77e`, preceded by `5ec2286`)
- Python 3.10.12 and pyprecice v3.0.0.0
- OpenFOAM v2312 and OpenFOAM-preCICE v1.3.0
- SU2 v7.5.1 and SU2-preCICE from commit `ac01cac` (effectively `91c3505`)
- CalculiX v2.20 and CalculiX-preCICE v2.20.1
- deal.II 9.4.0 and deal.II-preCICE from commit `2ab217d`
- FEniCS 2019.2.0.13.dev0 and fenicsprecice v2.1.0
- Nutils 8.6
- DUNE 2.8 and DUNE-preCICE from commit `e9fa630`
- solids4foam v2.0 on OpenFOAM v2012

Running on a different system is expected to have differences in the last digits.
In order to use these results as part of regression tests, choose a higher tolerance.
