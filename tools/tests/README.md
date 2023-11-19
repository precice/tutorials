---
title: preCICE system tests
permalink: dev-docs-system-tests.html
sidebar: docs_sidebar
keywords: pages, development, tests
summary: "Test complete simulations combining preCICE components of specific versions."
---

The tutorials repository hosts cases that need multiple components from the preCICE ecosystem to run. This directory provides tools that can automatically run complete simulations, using different versions of each component, and compare the results to references. While the main purpose is to run complete tests in the continuous integration workflows of preCICE, you can also run these tests on your laptop.

## Running the system tests

**! Warning: still under development !**

The main workflow for the user is executing the `systemtests.py` script. Depending on the options given to the script, it reads in the respective metadata files and generates `docker-compose.yaml` files that can start a fully-defined coupled simulation.

### Running specific test suites

To test a certain test-suite defined in `tests.yaml`, use:

```bash
python3 systemtests.py --suites=openfoam-adapter-release,<someothersuite>
```

To discover all tests, use `python print_test_suites.py`.
To be able to fill in the right case tuple into the `tests.yaml`, you can use the `python3 print_case_combinations.py` script.

## Adding new tests

### Adding tutorials

In order for the systemtests to pick up the tutorial we need to define a `metadata.yaml` in the folder of the tutorial. There are a few `metadata.yaml` already present to get inspiration from. You can also have a look at the implementation details but normally the currently available ones should be easy to adopt. You can check your metadata parsing by `python print_metadata.py` and `python print_case_combinations.py`

### Adding Testsuites

To add a testsuite just open the `tests.yaml` file and use the output of `python print_case_combinations.py` to add the right case combinations you want to test. Note that you can specify a `reference_result` which is not yet present. The `generate_reference_data.py` will pick that up and create it for you.
Note that its important to carefully check the paths of the `reference_result` in order to not have typos in there. Also note that same cases in different testsuites should use the same `reference_result`.

### Generate reference results

Since we need data to compare against, you need to run `python generate_reference_data.py`. This process might take a while.
Please include the generated reference results in the pull request as they are strongly connected to the new testsuites.

## Implementation details

Each tutorial contains automation scripts (mainly `run.sh` and `clean.sh`), as well as metadata (`metadata.yaml`). The metadata file describes the available cases, how to run them, as well as their dependencies. A central `tests.yaml` file in this directory defines test suites, which execute different combinations of cases. The Python script `systemtests.py` executes the tests, allowing to filter for specific components or test suites.

Let's dive deeper into some of these aspects.

### Metadata

Every tutorial contains a file called `metadata.yaml` describing some important properties of the tutorial. For example:

```yaml
name: Elastic tube 3D
path: elastic-tube-3d
url: https://precice.org/tutorials-elastic-tube-3d.html

participants:
  - Fluid 
  - Solid 

cases:
  fluid-openfoam:
    participant: Fluid
    directory: ./fluid-openfoam
    run: ./run.sh 
    component: openfoam-adapter
  
  solid-calculix:
    participant: Solid
    directory: ./solid-calculix
    run: ./run.sh 
    component: calculix-adapter 
  
  solid-fenics:
    participant: Solid
    directory: ./solid-fenics
    run: ./run.sh
    component: fenics-adapter 
```

Description:

- `name`: A human-readable, descriptive name
- `path`: Where the tutorial is located, relative to the tutorials repository
- `url`: A web page with more information on the tutorial
- `participants`: A list of preCICE participants, typically corresponing to different domains of the simulation
- `cases`: A list of solver configuration directories. Each element of the list includes:
  - `participant`: Which participant this solver case can serve as
  - `directory`: Where the case directory is located, relative to the tutorial directory
  - `run`: Command that executes the tutorial
  - `component`: Component or list of components that this case depends upon (typically an adapter)

### Components

The components mentioned in the Metadata are defined in the central `components.yaml` file. This file also specifies some arguments and their default values. These arguments can be anything, but most often they are version-related information. For example, the version of the OpenFOAM library used by the openfoam-adapter component, or the openfoam-adapter component version itself. For example:

```yaml
openfoam-adapter:
  repository: https://github.com/precice/openfoam-adapter
  template: component-templates/openfoam-adapter.yaml
  build_arguments: # these things mean something to the docker-service
    OPENFOAM_EXECUTABLE:
      options: ["openfoam2112"]
      description: exectuable of openfoam to use
      default: "openfoam2112"
    PRECICE_TAG:
      description: Version of preCICE to use
      default: "latest"
    OPENFOAM_ADAPTER_REF:
      description: Reference/tag of the actual OpenFOAM adapter
      default: "master"
```

This `openfoam-adapter` component has the following attributes:

- `repository`: URL to the Git project (without the `.git` extension), used to fetch the component
- `template`: A template for a Docker Compose service of this component
- `build_arguments`: Arguments passed to the Docker Compose service (arbitrary)

#### Naming schema for build_arguments

Since the docker containers are still a bit mixed in terms of capabilities and support for different build_argument combinations the following rules apply:

- A build_argument ending in **_TAG** means that an image of some kind needs to be available with that tag.
- A build_argument ending in **_REF** means that it refers to a git reference (like a branch or commit) beeing used to build the image.
- All other build_arguments are free of rules and up to the container maintainer.

### Component templates

Templates for defining a Docker Compose service for each component are available in `component-templates/`. For example:

```yaml
image: precice/fenics-adapter:{{ build_arguments["FENICS_ADAPTER_REF"] }}
depends_on:    
  prepare:
    condition: service_completed_successfully
volumes:
  - {{ run_directory }}:/runs
command: >
  /bin/bash -c "id && 
  cd '/runs/{{ tutorial_folder }}/{{ case_folder }}' &&
  {{ run }} | tee {{ case_folder }}.log 2>&1"
```

This template defines:

- `image`: The base Docker image for this component, including a Git reference (tag), provided to the template as argument (e.g., by the `systemtests.py` script).
- `depends_on`: Other services this service depends upon, typically a preparation service that fetches all components and tutorials.
- `volumes`: Directories mapped between the host and the container. Apart from directories relating to the users and groups, this also defines where to run the cases.
- `command`: How to run a case depending on this component, including how and where to redirect any screen output.

### Tests

Concrete tests are specified centrally in the file `tests.yaml`. For example:

```yaml
test_suites:
  openfoam_adapter_pr:
    tutorials:
      - path: flow-over-heated-plate
        case_combination:
          - fluid-openfoam
          - solid-openfoam
        reference_result: ./flow-over-heated-plate/reference-data/fluid-openfoam_solid-openfoam.tar.gz
  openfoam_adapter_release:
    tutorials:
      - path: flow-over-heated-plate
        case_combination:
          - fluid-openfoam
          - solid-openfoam
        reference_result: ./flow-over-heated-plate/reference-data/fluid-openfoam_solid-openfoam.tar.gz
      - path: flow-over-heated-plate
        case_combination:
          - fluid-openfoam
          - solid-fenics
        reference_result: ./flow-over-heated-plate/reference-data/fluid-openfoam_solid-fenics.tar.gz
```

This defines two test suites, namely `openfoam_adapter_pr` and `openfoam_adapter_release`. Each of them defines which case combinations of which tutorials to run.

### Generate Reference Results

In order to generate the reference results edit the `reference_versions.yaml` to match the required `build_arguments` otherwise passed via the cli.
Executing `generate_reference_data.py` will then generate a the following files:

- all distinct `.tar.gz` defined in the `tests.yaml`
- a `reference_results.md` in the tutorial folder describing the arguments used and a sha-1 hash of the `tar.gz` archive.

The reference result archive will later be unpacked again during the systemtest and compared using `fieldcompare`
Please note that these files should always be kept in the git lfs.
