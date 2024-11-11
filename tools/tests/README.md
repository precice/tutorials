---
title: preCICE system tests
permalink: dev-docs-system-tests.html
sidebar: docs_sidebar
keywords: pages, development, tests
summary: "Test complete simulations combining preCICE components of specific versions."
---

The tutorials repository hosts cases that need multiple components from the preCICE ecosystem to run. This directory provides tools that can automatically run complete simulations, using different versions of each component, and compare the results to references. While the main purpose is to run complete tests in the continuous integration workflows of preCICE, you can also run these tests on your laptop.

## Running the system tests

The main workflow for the user is executing the `systemtests.py` script. Depending on the options given to the script, it reads in the respective metadata files and generates `docker-compose.yaml` files that can start a fully-defined coupled simulation.

### Running the tests for a preCICE release

Workflow for the preCICE v3 release testing:

1. Collect the Git commits/tags of all components you want to test. The caching mechanism cannot detect changes based on branch names. The same effect might be encountered when rebasing and force-pushing the release branch.
2. In your terminal, navigate to the tutorials repository
3. Trigger the GitHub Actions Workflow. Until we merge the workflow to develop, this can only happen via the [GitHub CLI](https://cli.github.com/):

    ```bash
    gh workflow run run_testsuite_manual.yml -f suites=release_test -f build_args="PRECICE_REF:v3.1.1,PRECICE_PRESET:production-audit,OPENFOAM_ADAPTER_REF:v1.3.0,PYTHON_BINDINGS_REF:v3.1.0,FENICS_ADAPTER_REF:v2.1.0,SU2_VERSION:7.5.1,SU2_ADAPTER_REF:64d4aff,TUTORIALS_REF:340b447" --ref=develop
    ```

4. Go to the tutorials [Actions](https://github.com/precice/tutorials/actions) page and find the running workflow
5. Check the status and the runtimes of each tutorial:

    - Very small build times mean that the test is using cached container layers
    - Most commonly, you will see tests failing with `Fieldcompare returned non zero exit code`. You will need to check the logs, but if the fieldcompare time is significant, this typically means that the numerical results differ above the tolerance (the test works!).

6. Download the build artifacts from Summary > runs.

    - In there, you may want to check the `stdout.log` and `stderr.log` files.
    - The produced results are in `precice-exports/`, the reference results in `reference-results-unpacked`.
    - Compare using, e.g., ParaView or [fieldcompare](https://gitlab.com/dglaeser/fieldcompare): `fieldcompare dir precice-exports/ reference/`. The `--diff` option will give you `precice-exports/diff_*.vtu` files, while you can also try different tolerances with `-rtol` and `-atol`.

### Running specific test suites

To test a certain test-suite defined in `tests.yaml`, use:

```bash
python3 systemtests.py --suites=fenics_test,<someothersuite>
```

To discover all tests, use `python print_test_suites.py`.

To be able to fill in the right case tuple into the `tests.yaml`, you can use the `python3 print_case_combinations.py` script.

## Running the system tests on GitHub Actions

Go to Actions > [Run Testsuite (manual)](https://github.com/precice/tutorials/actions/workflows/run_testsuite_manual.yml) to see this workflow.

After bringing these changes to `master`, the manual triggering option should be visible on the top right. Until that happens, we can only trigger this workflow manually from the [GitHub CLI](https://github.blog/changelog/2021-04-15-github-cli-1-9-enables-you-to-work-with-github-actions-from-your-terminal/):

```shell
gh workflow run run_testsuite_manual.yml -f suites=fenics_test --ref=develop
```

Another example, to use the latest releases and enable debug information of the tests:

```shell
gh workflow run run_testsuite_manual.yml -f suites=fenics_test -f build_args="PRECICE_REF:v3.1.1,PRECICE_PRESET:production-audit,OPENFOAM_ADAPTER_REF:v1.3.0,PYTHON_BINDINGS_REF:v3.1.0,FENICS_ADAPTER_REF:v2.1.0,SU2_VERSION:7.5.1,SU2_ADAPTER_REF:64d4aff,TUTORIALS_REF:340b447" -f loglevel=DEBUG --ref=develop
```

where the `*_REF` should be a specific [commit-ish](https://git-scm.com/docs/gitglossary#Documentation/gitglossary.txt-aiddefcommit-ishacommit-ishalsocommittish).

Example output:

```text
Run cd tools/tests
  cd tools/tests
  python systemtests.py --build_args=PRECICE_REF:v3.1.1,PRECICE_PRESET:production-audit,OPENFOAM_ADAPTER_REF:v1.3.0,PYTHON_BINDINGS_REF:v3.1.0,FENICS_ADAPTER_REF:v2.1.0 --suites=fenics_test --log-level=DEBUG
  cd ../../
  shell: /usr/bin/bash -e {0}
INFO: About to run the following systemtest in the directory /home/precice/runners_root/actions-runner-tutorial/_work/tutorials/tutorials/runs:
 [Flow over heated plate (fluid-openfoam, solid-fenics)]
INFO: Started running Flow over heated plate (fluid-openfoam, solid-fenics),  0/1
DEBUG: Checking out tutorials master before copying
From https://github.com/precice/tutorials
 * [new branch]      master     -> master
DEBUG: Building docker image for Flow over heated plate (fluid-openfoam, solid-fenics)
DEBUG: Running tutorial Flow over heated plate (fluid-openfoam, solid-fenics)
DEBUG: Running fieldcompare for Flow over heated plate (fluid-openfoam, solid-fenics)
DEBUG: extracting /home/precice/runners_root/actions-runner-tutorial/_work/tutorials/tutorials/flow-over-heated-plate/reference-results/fluid-openfoam_solid-fenics.tar.gz into /home/precice/runners_root/actions-runner-tutorial/_work/tutorials/tutorials/runs/flow-over-heated-plate_fluid-openfoam-solid-fenics_2023-11-19-211723/reference_results
Using log-level: DEBUG
+---------------------------------------------------------+---------+-------------------+-----------------+-----------------------+
| systemtest                                              | success | building time [s] | solver time [s] | fieldcompare time [s] |
CRITICAL: Fieldcompare returned non zero exit code, therefore Flow over heated plate (fluid-openfoam, solid-fenics) failed
INFO: Running Flow over heated plate (fluid-openfoam, solid-fenics) took 280.5861554039875 seconds
ERROR: Failed to run Flow over heated plate (fluid-openfoam, solid-fenics)
+---------------------------------------------------------+---------+-------------------+-----------------+-----------------------+
| Flow over heated plate (fluid-openfoam, solid-fenics)   |    0    |      271.80       |      5.60       |         2.42          |
+---------------------------------------------------------+---------+-------------------+-----------------+-----------------------+
```

In this case, building and running seems to work out, but the tests fail because the results differ from the reference results. This may be incorrect, as the previous step may have silently failed.

## Understanding what went wrong

The easiest way to debug a systemtest run is first to have a look at the output written into the action on GitHub.
If this does not provide enough hints, the next step is to download the generated `runs` artifact. Note that by default this will only be generated if the systemtests fail.
Inside the archive, a test-specific subfolder like `flow-over-heated-plate_fluid-openfoam-solid-fenics_2023-11-19-211723` contains two log files: a `stderr.log` and `stdout.log`. This can be a starting point for a further investigation.

## Adding new tests

### Adding tutorials

In order for the systemtests to pick up the tutorial we need to define a `metadata.yaml` in the folder of the tutorial. There are a few `metadata.yaml` already present to get inspiration from. You can also have a look at the implementation details but normally the currently available ones should be easy to adopt. You can check your metadata parsing by `python print_metadata.py` and `python print_case_combinations.py`

### Adding testsuites

To add a testsuite just open the `tests.yaml` file and use the output of `python print_case_combinations.py` to add the right case combinations you want to test. Note that you can specify a `reference_result` which is not yet present. The `generate_reference_data.py` will pick that up and create it for you.
Note that its important to carefully check the paths of the `reference_result` in order to not have typos in there. Also note that same cases in different testsuites should use the same `reference_result`.

### Generate reference results

Since we need data to compare against, you need to run `python generate_reference_data.py`. This process might take a while.
Please include the generated reference results in the pull request as they are strongly connected to the new testsuites.

## Adding repositories

If you want to trigger a testsuite from a new repository, you need to add a workflow file to that repository (under `.github/workflows/`). You can, for example, copy and adjust [the one from the OpenFOAM adapter](https://github.com/precice/openfoam-adapter/blob/develop/.github/workflows/system-tests.yaml). Then, you need a new label to trigger the workflow (e.g. `Issues`->`Labels`->`New label`). Last, in the [preCICE organization settings](https://github.com/organizations/precice/settings), you need to add the new repository to the action secret `WORKFLOW_DISPATCH_TOKEN` and to the default actions runner group. Adding the new label directly to the pull request in which you add the workflow should already trigger the system tests, compare the [pull request in the `precice` repository](https://github.com/precice/precice/pull/2052).

## Implementation details

Each tutorial contains automation scripts (mainly `run.sh` and `clean.sh`), as well as metadata (`metadata.yaml`). The metadata file describes the available cases, how to run them, as well as their dependencies. A central `tests.yaml` file in this directory defines test suites, which execute different combinations of cases. The Python script `systemtests.py` executes the tests, allowing to filter for specific components or test suites.

Let's dive deeper into some of these aspects.

### General architecture

Each tutorial directory contains a metadata file, describing which participants each case directory implements, and how to run it.

A list of tests describes all tests to be executed, grouped by test suites. Each test is a combination of tutorial cases.

Test steps include modifying the tutorial configuration files for the test system, building the Docker containers used by the respective Docker Compose service of each component, and comparing results to reference results using fieldcompare.

Tests are executed by the `systemtests.py` script, which starts the Docker Compose. This can be executed locally, and it is the same script that GitHub Actions also execute.

The multi-stage Docker build allows building each component separately from the same Dockerfile, while Docker reuses cached layers. The Docker Compose services consider GitHub Actions Cache when building the services, although the cache is currently only updated, but not hit (see https://github.com/precice/tutorials/pull/372#issuecomment-1748335750).

### File structure

Metadata and workflow/script files:

- `.github/workflows/`
  - `run_testsuite_workflow.yml`: workflow for running the tests, triggered by other workflows (e.g., other repositories)
  - `run_testsuite_manual.yml`: manual triggering front-end for `run_testsuite_workflow.yml`
- `flow-over-a-heated-plate/`
  - `fluid-openfoam/`
    - `run.sh`: describes how to execute the respective case
  - `solid-fenics/`
  - `solid-openfoam/`
  - ...
  - `metadata.yml`: describes each case directory (which participant, which component, which script to run, ...)
- `tools/tests/`
  - `component-templates/`: jinja2 templates for Docker Compose services for the components
    - `calculix-adapter.yaml`
    - `fenics-adapter.yaml`
    - `openfoam-adapter.yaml`
    - ...
  - `dockerfiles/ubuntu_2204/`
    - Dockerfile: a multi-stage build Dockerfile that defines how to build each component, in a layered approach
  - `docker-compose.template.yaml`: Describes how to prepare each test (Docker Componse service template)
  - `docker-compose.field_compare.template.yaml`: Describes how to compare results with fieldcompare (Docker Compose service template)
  - `components.yaml`: Declares the available components and their parameters/options
  - `reference_results.metadata.template`: Template for reporting the versions used to generate the reference results
  - `reference_versions.yaml`: List of arguments to use for generating the reference results
  - `tests.yaml`: Declares the available tests, grouped in test suites

User-facing tools:

- `tools/tests/`
  - `systemtests.py`: Executes the system tests, starting Docker Compose services of each required component (after building them), running each test, and comparing the results to reference results.
  - `print_test_suites.py`: Prints the available tests.
  - `print_metadata.py`: Prints the metadata of each tutorial that contains a `metadata.yaml` file.
  - `print_case_combinations.py`: Prints all possible combinations of tutorial cases, using the `metadata.yaml` files.
  - `build_docker_images.py`: Build the Docker images for each test
  - `generate_reference_results.py`: Executes the system tests with the versions defined in `reference_versions.yaml` and generates the reference data archives, with the names described in `tests.yaml`. (should only be used by the CI Pipeline)

Implementation scripts:

- `tools/tests/`
  - `systemtests.py`: Main entry point
  - `requirements.txt`: Dependencies (jinja2, pyyaml)
  - `metadata_parser/`: Reads the YAML files into Python objects (defines the schema)
  - `systemtests/`: Main implementation classes
    - `Systemtest.py`
    - `SystemtestArguments.py`
    - `TestSuite.py`

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
  build_arguments:
    PRECICE_REF:
      description: Version of preCICE to use
      default: "main"
    PRECICE_PRESET:
      description: CMake preset of preCICE
      default: "production-audit"
    PLATFORM:
      description: Dockerfile platform used
      default: "ubuntu_2204"
    TUTORIALS_REF:
      description: Tutorial git reference to use
      default: "master"
    OPENFOAM_EXECUTABLE:
      options: ["openfoam2306","openfoam2212","openfoam2112"]
      description: exectuable of openfoam to use
      default: "openfoam2306"
    OPENFOAM_ADAPTER_REF:
      description: Reference/tag of the actual OpenFOAM adapter
      default: "master"
```

This `openfoam-adapter` component has the following attributes:

- `repository`: URL to the Git projects
- `template`: A template for a Docker Compose service of this component
- `build_arguments`: Arguments passed to the Docker Compose service (arbitrary)

#### Naming schema for build_arguments

Since the docker containers are still a bit mixed in terms of capabilities and support for different build_argument combinations the following rules apply:

- A build_argument ending in **_REF** means that it refers to a git commit-ish (like a tag or commit) beeing used to build the image. Its important to not use branch names here as we heavily rely on dockers build cache to speedup things. But since the input variable to the docker builder will not change, we might have wrong cache hits.
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
  {{ run }} | tee system-tests_{{ case_folder }}.log 2>&1"
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
        reference_result: ./flow-over-heated-plate/reference-results/fluid-openfoam_solid-openfoam.tar.gz
  openfoam_adapter_release:
    tutorials:
      - path: flow-over-heated-plate
        case_combination:
          - fluid-openfoam
          - solid-openfoam
        reference_result: ./flow-over-heated-plate/reference-results/fluid-openfoam_solid-openfoam.tar.gz
      - path: flow-over-heated-plate
        case_combination:
          - fluid-openfoam
          - solid-fenics
        reference_result: ./flow-over-heated-plate/reference-results/fluid-openfoam_solid-fenics.tar.gz
```

This defines two test suites, namely `openfoam_adapter_pr` and `openfoam_adapter_release`. Each of them defines which case combinations of which tutorials to run.

### Generate Reference Results

#### via GitHub workflow (recommended)

The preferred way of adding reference results is via the manual triggerable `Generate reference results (manual)` workflow. This takes two inputs:

- `from_ref`: branch where the new test configuration (e.g added tests, new reference_versions.yaml) is
- `commit_msg`: commit message for adding the reference results into the branch

The workflow will checkout the `from_ref`, take the status of the systemtests of that branch and execute `python generate_reference_results.py`, upload the LFS objects into the self-hosted LFS server and add a commit with `commit_msg` onto the `from_ref` branch.

#### manually

In order to generate the reference results edit the `reference_versions.yaml` to match the required `build_arguments` otherwise passed via the cli.
Executing `generate_reference_results.py` will then generate the following files:

- all distinct `.tar.gz` defined in the `tests.yaml`
- a `reference_results.md` in the tutorial folder describing the arguments used and a sha-1 hash of the `tar.gz` archive.

The reference result archive will later be unpacked again during the systemtest and compared using `fieldcompare`
Please note that these files should always be kept in the git lfs.
