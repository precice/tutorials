---
title: Partitioned Pipe Multiscale
permalink: tutorials-partitioned-pipe-multiscale.html
keywords: OpenFOAM, python
summary: The 1D-3D Partitioned Pipe is a simple geometric multiscale case, coupling two pipes with different dimensions. 
---

{% note %}
Get the [case files of this tutorial](https://github.com/precice/tutorials/tree/master/partitioned-pipe-multiscale). Read how in the [tutorials introduction](https://www.precice.org/tutorials.html).
{% endnote %}

## Setup

We exchange velocity data from the upstream 1D to the downstream 3D participant and for the pressure data vice versa. The config looks as follows:

![Config visualization](images/tutorials-partitioned-pipe-multiscale-config.png)

## How to run

In two different terminals execute

```bash
cd fluid1d-python && ./run.sh
```

```bash
cd fluid3d-openfoam && ./run.sh
```

## Results

Visualizing the results in ParaView, we see an established laminar profile at the inlet of the 3D participant.

![Expected result](images/tutorials-partitioned-pipe-multiscale-profile.png)
