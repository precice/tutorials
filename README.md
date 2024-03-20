# preCICE tutorials

> [!IMPORTANT]  
> This repository is aimed for development purposes and the default branch is `develop`. If you want to use the tutorials, switch to the [`master` branch](https://github.com/precice/tutorials/tree/master) or download the latest [release](https://github.com/precice/tutorials/releases).

This repository contains ready-to-run tutorial cases for the coupling library [preCICE](https://precice.org/).
The purpose of these cases is not to teach you how to use preCICE from scratch, but to serve as starting points for setting up similar simulation cases, as well as test cases. Read more on our [preCICE tutorials](https://precice.org/tutorials.html) documentation section.

As a general rule, you can start each participant from inside their `<tutorial>/<participant>-<solver>` using `./run.sh`. Look into these short scripts and copy the parts you need for your new case. Before running again, execute the cleaning scripts you can find at each level, to clean from this point and deeper.

Contributions to this repository are very welcome. Please refer to the page [Contribute to preCICE](https://precice.org/community-contribute-to-precice.html) for a few guidelines and hints to help you in this direction.

Note that we use [Git LFS](https://git-lfs.com/) to version reference results. These will appear as seemingly empty files containing URLs if you don't have Git LFS installed (optional, mainly useful for our system tests).
