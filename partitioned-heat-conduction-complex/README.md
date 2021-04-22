---
title: Partitioned heat conduction (complex setup)
permalink: tutorials-partitioned-heat-conduction-complex.html
keywords: FEniCS, Heat conduction
summary: This tutorial is the advanced version of the "partitioned heat conduction" tutorial, showcasing more advanced features and geometries.
---

{% include note.html content="Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-complex). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html)." %}

## Setup

This case is an advanced version of `partitioned-heat-conduction`. Some advanced features offered by this case:

* Geometries may be chosen arbitrarily. One possibility is to use a circle and a rectangular plate with a hole, but you can also provide your own geometry, if you want.
* You may combine arbitrary mesh resolutions at the coupling interface.
* Nearest projection mapping is used.
* The Dirichlet and Neumann participants may be swapped arbitrarily.
* The exchanged temperature is still scalar valued, but the heat flux is vector valued.
* You can decide to use a time dependent heat flux and right-hand side to make the problem more challenging.

## Available solvers and dependencies

See `partitioned-heat-conduction`, only `fenics` is provided as a solver.

## Running the simulation

See `partitioned-heat-conduction`. The additional featured mentioned above can be activated via command line arguments. Please run `python3 fenics/heat.py --help` for a full list of provided arguments.
