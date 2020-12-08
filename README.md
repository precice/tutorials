# preCICE Tutorials

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](http://www.precice.org/).

You may find step-by-step instructions for each case in the [preCICE wiki](https://github.com/precice/precice/wiki). *More tutorials come with each adapter* and you can also find them in the wiki.


Each folder in this directory represents a stand-alone tutorial containing:
* instructions in the `docs` folder.
* multiple participant folders with the format `<participant>_<solver>`. Pick a solver for each participant and start them to run the simulation.


The website import relies on the following rules:
* tutorials use dashes to separate words: `flow-over-plate`
* a tutorial contains a `docs` folder.
* `docs/instructions.md` contains a Jekyll front-matter and a `permalink: /tutorials-flow-over-plate.html`.
* `docs/images/` is optional and all content will be available in `images/tutorials/flow-over-plate/`. (Use this prefix to include images in `instructions.md`)
