# Systemtests
## What do they solve?
Since the last Systemtest got deprecated at some point a new precice release involves a lot of manual system testing to ensure stable operations between the different adapters. 

## Why is this in the tutorial repository?
Most of the cases we currently test new precice releases on are the tutorial cases. They cover quite a broad spectrum of different solvers and domains. Putting everything into a seperate repository would make maintance harder as you would split up things that belong together. 

## How do the systemtests work
Since we rely heavily on the tutorials to test the preCICE on a system level, we added tutorial metadata to descibe the tutorials in a systemmatic way. This enables us to programatically run them to test the software. 
In the following we will roughly introduce some bigger picture concepts used in the metadata and needed for the systemtests. 

### Tutorial description
Every tutorial contains a file called `metadata.yaml` describing some important properties of the tutorial. In the following most of the possible entries into this `yaml` file will be explained.
#### Participants
Participants are the basic builing blocks of a coupled simulation. They describe what will be coupled together. One very common participant couple would be fluid participant and a solid participant doing an fluid-structure interaction simulation.
#### Components
Additionally the `metadata.yaml` contains a list of available **components** that will be responsible for the actual solution of the isolated problem. An example component would be the openfoam-adapter or the fenics-adapter. These components are also described globally in a `components.yaml` file. In that file it also specifies some **arguments**, and their defaults given to the **components**. These arguments can be anything, but most often they are version related information. For example the version of the OpenFOAM library used by the openfoam-adapter or the openfoam-adapter version itself. 
#### Cases
Lastly we have a case. A Case is the instantiaton of a **participant** using the code of a **component** to do so. Cases correspond to the different subfolders in the tutorial directory. Therefore you need to specify the folder the case is contained and provide it with a `run` command. 


## How do i execute them?
! Warning still a WIP

Therefore the workflow for the user is to just execute the `systemtests.py` file. Depending on the options given to the file i will then read in all the metada files and generate the appropriate `docker-compose.yaml` files. 

### By Test-Suite
To test a certain test-suite defined in `tests.yaml` just use `python systemtests.py --suites=openfoam-adapter-release,<someothersuite>`. 
To discover all tests use `python print_test_suites.py`. 
To be able to fill in the right case tuple into the `tests.yaml` you can use the `python print_case_combinations.py` script. 


### By component
! Do not use, this will be deprecated at some point
To test the current state, which only supports openfoam run:
` python systemtests.py --components=openfoam-adapter --params=openfoam-version=v2012,`

