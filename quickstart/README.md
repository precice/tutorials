---
title: Quickstart
permalink: quickstart.html
keywords: tutorial, quickstart
summary: "Install preCICE on Linux (e.g. via a Debian package), and then couple two OpenFOAM solvers with the OpenFOAM-preCICE adapter."
layout: "page"
comments: false
search: true
sidebar: nil
topnav: topnav
toc: false
---



## Start here

1. To get a feeling what preCICE does, watch a [short presentation](https://www.youtube.com/watch?v=FCv2FNUvKA8), a [longer training session](https://www.youtube.com/watch?v=FCv2FNUvKA8), or [click through a tutorial in your browser](http://run.precice.org/).
2. Get and install preCICE. For Linux systems, this is pretty easy. Just pick what suits you best on [this overview page](installation-overview.html). Facing any problems? [Ask for help](community-channels.html).
    - For example, [download](https://github.com/precice/precice/releases/latest) and install our binary package for Ubuntu 20.04 (Focal Fossa) by clicking on it or using the following commands:
    ```shell
    wget https://github.com/precice/precice/releases/download/v2.1.1/libprecice2_2.1.1_focal.deb
    sudo apt install ./libprecice2_2.1.1_focal.deb
    ```
3. Build and couple two [C++ solverdummies](https://github.com/precice/precice/tree/master/examples/solverdummies/cpp).
4. You probably want to couple a solver you are already using, such as OpenFOAM. Since many of our tutorials use it, [install OpenFOAM](adapter-openfoam-support.html).
5. Download and install the [OpenFOAM-preCICE adapter](adapter-openfoam-get.html).
6. [Couple OpenFOAM with OpenFOAM](https://github.com/precice/openfoam-adapter/wiki/Tutorial-for-CHT:-Flow-over-a-heated-plate). This testcase is part of the OpenFOAM adapter.

## What's next?

To become a preCICE pro:

* Get an overview of the [preCICE docs](docs.html).
* See what users talk about in the [preCICE forum](https://precice.discourse.group/).
* Run [tutorials with other coupled solvers](https://github.com/precice/precice/wiki#2-getting-started---tutorials).
* Watch some [preCICE videos](https://www.youtube.com/channel/UCxZdSQdmDrheEqxq8g48t6A).
* Register to our [virtual preCICE Workshop 2021](precice-workshop-2021.html).
* Find out how to [couple your own solver](couple-your-code-prerequisites.html).
* Tell us [your story](community-projects.html).


{% include important.html content="We are currently working on a really nice getting started page. Till then, we hope these tips already help." %}
