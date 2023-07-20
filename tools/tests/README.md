# preCICE system tests

The tutorials repository hosts cases that need multiple components from the preCICE ecosystem to run. This directory provides tools that can automatically run complete simulations, using different versions of each component, and compare the results to references. While the main purpose is to run complete tests in the continuous integration workflows of preCICE, you can also run these tests on your laptop.

## Running the system tests

**! Warning: still under development !**

The main workflow for the user is executing the `systemtests.py` script. Depending on the options given to the script, it reads in the respective metadata files and generates `docker-compose.yaml` files that can start a fully-defined coupled simulation.

### Running specific test suites

To test a certain test-suite defined in `tests.yaml`, use:

```bash
python systemtests.py --suites=openfoam-adapter-release,<someothersuite>
```

To discover all tests, use `python print_test_suites.py`.
To be able to fill in the right case tuple into the `tests.yaml`, you can use the `python print_case_combinations.py` script.

### Testing specific components

**! Do not use, this will be deprecated at some point !**

To test the current state, which only supports OpenFOAM, run:

```bash
python systemtests.py --components=openfoam-adapter --params=openfoam-version=v2012
```

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
  build-arguments: # these things mean something to the docker-service
    OPENFOAM_EXECUTABLE:
      options: ["openfoam2306"]
      description: exectuable of OpenFOAM to use
      default: "openfoam2306"
    PRECICE_TAG:
      description: Version of precice to use
      default: "latest"
    OPENFOAM_ADAPTER_TAG:
      description: Ref of the actual OpenFOAM adapter
      default: "latest"
```

This `openfoam-adapter` component has the following attributes:

- `repository`: URL to the Git project (without the `.git` extension), used to fetch the component
- `template`: A template for a Docker Compose service of this component
- `build-arguments`: Arguments passed to the Docker Compose service (arbitrary)

### Component templates

Templates for defining a Docker Compose service for each component are available in `component-templates/`. For example:

```yaml
image: "ghcr.io/precice/openfoam-adapter:{{ params["openfoam-adapter-ref"] }}"
user: ${MY_UID}:${MY_GID}
depends_on:    
  prepare:
    condition: service_completed_successfully
volumes:
  - /etc/passwd:/etc/passwd:ro
  - /etc/group:/etc/group:ro
  - {{ run_directory }}:/runs
command: >
  /bin/bash -c "id && 
  cd '/runs/{{ tutorial_folder }}/{{ case_folder }}' &&
  {{ params["openfoam-exectuable"] }} {{ run }} | tee {{ case_folder }}.log 2>&1"
```

This template defines:

- `image`: The base Docker image for this component, including a Git reference (tag), provided to the template as argument (e.g., by the `systemtests.py` script).
- `user`: The user executing the service, provided to the template as argument (e.g., by the `systemtests.py` script). This is important to ensure that the result files created are owned by the user that started the system tests.
- `depends_on`: Other services this service depends upon, typically a preparation service that fetches all components and tutorials.
- `volumes`: Directories mapped between the host and the container. Apart from directories relating to the users and groups, this also defines where to run the cases.
- `command`: How to run a case depending on this component, including how and where to redirect any screen output.

### Tests

Concrete tests are specified centrally in the file `tests.yaml`. For example:

```yaml
test-suites:
  openfoam-adapter-pr:
    tutorials:
      flow-over-heated-plate:
        cases:
          - (fluid-openfoam, solid-openfoam)
          - (fluid-openfoam, solid-nutils)
      elastic-tube-3d:
        cases:
          - (fluid-openfoam, solid-fenics)
  openfoam-adapter-release:
    tutorials:
      flow-over-heated-plate:
        cases:
          - (fluid-openfoam, solid-openfoam)
```

This defines two test suites, namely `openfoam-adapter-pr` and `openfoam-adapter-release`. Each of them defines which case combinations of which tutorials to run.
