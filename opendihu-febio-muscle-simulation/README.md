# Bachelor-Forschungsprojekt: Comparison and coupling of the neuromuscular simulation software OpenDiHu with the biomechanics solver FEBio

This is the repository for the Bachelor-Forschungsprojekt of Paul Arlt, Luis Morgenstern, Silas Natterer and Jan Stein. It encompasses all the code and documentation related to our research project on modeling skeletal muscles using the FEBio and OpenDiHu simulation frameworks. For a more detailed description of the project and its results, we refer to [our report](report.pdf).

We've organized the repository into several folders, each representing a distinct aspect of our project:

- [bfp_plugin](bfp_plugin): In this folder, you'll find the code for our custom FEBio plugin, originally forked from <https://github.com/precice/febio-adapter>. This adapter enables the coupling of FEBio and OpenDiHu through the library preCICE.

- [comparison_opendihu_febio](comparison_opendihu_febio): This folder encompasses the simulation cases associated with the comparative analysis of the mechanical simulation capabilities between OpenDiHu and FEBio.

- [active_contraction_cases](active_contraction_cases): This folder is dedicated to the simulation of active fiber contractions in OpenDiHu with and without preCICE. This also contains the main result of our project: A hybrid simulation, combining the mechanical solver of FEBio with the fast monodomain solver of OpenDiHu.
