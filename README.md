# preCICE Tutorials

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](http://www.precice.org/).
The purpose of these cases are not to teach you how to use preCICE, but to serve as starting points for setting up similar simulation cases.

Each folder in this directory represents a stand-alone tutorial case containing:
* a `README.md` with general instructions on how to run the simulation. This is the same content you see on the [preCICE tutorials](https://www.precice.org/tutorials.html) documentation section.
* multiple participant folders with the format `<participant>-<solver>`. Pick a solver for each participant and start them to run the simulation.
* a `clean.sh` script to clean each individual participant case, a `clean-tutorial.sh` to clean all participants of a tutorial case, a `clean-all.sh` to clean all tutorial cases.

As a general rule, you can start each participant from inside their `<participant>-<solver>` using `./run.sh`. Look into these short scripts and copy the parts you need for your new case.

## For contributors

The [preCICE tutorials](https://www.precice.org/tutorials.html) documentation section displays material from this repository. This website import relies on the following rules:
* tutorials use dashes to separate words `flow-over-heated-plate`.
* a tutorial contains a `REAMDE.md` with a Jekyll front-matter and a `permalink: /tutorials-flow-over-heated-plate.html`.
* a tutorial may contain a folder `images/` with images to be displayed in the `README.md`.  
  Images use the naming scheme `tutorials-flow-over-heated-plate-example.png`. Subdirectories in `images/` are not allowed.
* images includes in the `README.md` use relative links to `![example](images/tutorials-flow-over-heated-plate-example.png)`.
