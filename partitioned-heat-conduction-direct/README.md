---
title: Partitioned heat conduction (direct access setup)
permalink: tutorials-partitioned-heat-conduction-direct.html
keywords: Nutils, Heat conduction, Direct mesh access
summary: This tutorial is a modified version of the "partitioned heat conduction" tutorial showcasing direct mesh access.
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-direct). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

This case is a modified version of the [partitioned heat conduction tutorial](tutorials-partitioned-heat-conduction.html).

We use the [direct mesh access feature](couple-your-code-direct-access.html) to let the solvers compute the data mapping and not preCICE.

## Available solvers

Currently only `nutils` is provided as a solver.

## Running the simulation

Open two terminals and run:

```bash
cd nutils
./run.sh -d
```

and

```bash
cd fenics
./run.sh -n
```

See the [partitioned heat conduction tutorial](tutorials-partitioned-heat-conduction.html).
