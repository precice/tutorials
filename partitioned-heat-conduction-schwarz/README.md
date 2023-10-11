---
title: Partitioned heat conduction with Schwarz-type domain decomposition
permalink: tutorials-partitioned-heat-conduction-schwarz.html
keywords: FEniCS, Nutils, Heat conduction
summary: We solve a simple heat equation. The domain is partitioned and the coupling is established in an overlapping-Schwarz-type domain decomposition.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-schwarz). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Running the simulation

This tutorial is for FEniCS.

For choosing whether you want to run the left or right participant, please provide the following commandline input:

* `heat.py left` flag will run the left participant.
* `heat.py right` flag will run the right participant.
