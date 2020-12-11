# preCICE Tutorials

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](http://www.precice.org/).

You may find step-by-step instructions for each case in the [preCICE wiki](https://github.com/precice/precice/wiki). *More tutorials come with each adapter* and you can also find them in the wiki.


Each folder in this directory represents a stand-alone tutorial containing:
* a `README.md` with general instructions on how to run the simulation
* multiple participant folders with the format `<participant>_<solver>`. Pick a solver for each participant and start them to run the simulation.


The website import relies on the following rules:
* tutorials use dashes to separate words `flow-over-plate`.
* a tutorial contains a `REAMDE.md` with a Jekyll front-matter and a `permalink: /tutorials-flow-over-plate.html`.
* a tutorial may contain a folder `images/` with images to be displayed in the `README.md`.  
  Images use the naming scheme `tutorials-flow-over-plate-example.png`. Subdirectories in `images/` are not allowed.
* images includes in the `README.md` use relative links to `![example](images/tutorials-flow-over-plate-example.png)`.
