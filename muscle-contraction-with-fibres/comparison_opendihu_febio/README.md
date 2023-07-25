# Comparison OpenDiHu FEBio

This folder contains a common case for both softwares OpenDiHu and FEBio. In case of FEBio, the case is modeled using the SI unit system, as well as using the OpenDiHu unit system.

## Cases
The following cases are provided:
- [muscle_no_fibers](muscle_no_fibers): A muscle without any fibers, fixed on one end, pulled on the other end with a force increasing over time.
- [muscle_no_fibers_scaled](muscle_no_fibers_scaled): Present only for FEBio. Same as muscle_no_fibers but uses SI unit system to model the case.

## Running the OpenDiHu cases
To run an OpenDiHu case, you first have to compile it using
```bash
opendihu/scripts/shortcuts/sr.sh
```
After this you can navigate to *build_release* and run
```bash
./program_name ../settings.py
```

## Running the FEBio cases
To run a FEBio case, run
```bash
path_to_febio4_exe/febio4 case_name
```
from the *febio* folder.
