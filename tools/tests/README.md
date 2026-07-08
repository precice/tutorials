---
title: preCICE system tests
permalink: dev-docs-system-tests.html
aliases:
  - /dev-docs-system-tests.html
sidebar: docs_sidebar
keywords: pages, development, tests
summary: "Test complete simulations combining preCICE components of specific versions."
---

The tutorials repository hosts cases that need multiple components from the preCICE ecosystem to run. Te system tests automatically run tutorials combining the required components and compare the numerical results to references. While the main purpose is to run complete tests in the continuous integration workflows of preCICE, you can also run these tests on your laptop.

## Running

The main workflow for the user is executing the `systemtests.py` script, which is the same that the [GitHub Actions workflow](https://github.com/precice/tutorials/actions/workflows/system-tests-latest-components.yml) executes. Depending on the options given to the script, it reads in the respective metadata files and generates `docker-compose.yaml` files that can start a fully-defined coupled simulation. For arguments that are not provided, default values from `components.yaml` are used.

### Manual and nightly runs on GitHub

The [System tests (manual/nightly)](https://github.com/precice/tutorials/actions/workflows/system-tests-latest-components.yml) workflow executes the `release` test suite nightly and can also be triggered manually.

On the workflow page, click `Run workflow`. The default values will execute the `release` test suite using the latest `develop` branches of every component. If you want to override the version of some component, specify it in the respective field. Commit hashes, branches, and tags are all accepted. Branches and tags will get automatically resolved to their current commit on GitHub before starting any test, and all tests will use the same version of any common component.

The available test suites are found in [`tests.yaml`](https://github.com/precice/tutorials/blob/develop/tools/tests/tests.yaml) and common values are:

- `quickstart`, `elastic-tube-1d`, or any other tutorial (see [exceptions](https://github.com/precice/tutorials/issues/448))
- `openfoam-adapter`, `micro-manager`, `fmi-runner`, or similar test cases involving the respective component
- `precice` is a subset of cases that cover a range of preCICE features
- `release` is for all available but some very long or known to fail tests
- `extra` is for some longer tests
- `selected` and `system-tests-dev` for some special cases

The `Use workflow from` is a default option of GitHub Actions that concerns the GHA workflow file itself.

### Running from a pull request

Several repositories include a workflow that allows triggering the system tests by adding the `trigger-system-tests` label to the pull request. The event is triggered only at the moment of adding the label, so you need to remove the label and add it again if needed.

See the [system tests workflow in the OpenFOAM adapter](https://github.com/precice/openfoam-adapter/blob/develop/.github/workflows/system-tests.yaml) for an example.

### Running from the GitHub CLI

The [GitHub CLI](https://cli.github.com/) allows triggering workflows of a GitHub project. For example:

```bash
gh workflow run system-tests-latest-components.yml -f suites=release
```

More arguments are available, for example:

```bash
gh workflow run system-tests-latest-components.yml -f suites=release -f build_args="PLATFORM:ubuntu2404,PRECICE_REF:develop" -f log_level="DEBUG" --ref=develop
```

The `build_args` override the defaults set in `tools/tests/components.yaml`.

### Running locally

To run locally, you will need Docker, Docker Compose, and Python 3.

Navigate into the directory `tools/tests/` of the tutorials, make a Python virtual environment, and install the dependencies:

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements.txt
```

To test a certain test-suite defined in `tests.yaml`, use:

```bash
python systemtests.py --suites=quickstart
```

This will build and connect Docker containers and run tutorials in the `runs/` directory in the root of the tutorials.

To clean up at the end, you might want to run a `docker system prune -a` and remove the `runs/` directory to save space.

There are also some auxiliary scripts (run without arguments):

- `print_test_suites.py`: Print all test suites defined in `tests.yaml`
- `print_case_combinations.py`: Print all possible combinations of participants in a tutorial, using its `metadata.yaml`.

## Understanding the logs

For a local run, look into the `tutorials/runs/` directory, where you will find a time-stamped directory for every test executed.

Each of these directories includes the usual tutorial case files and logs, as well as:

1. `system-tests-build.log`: The logs of building the respective components.
2. `system-tests-run.log`: The logs of running the simulation (intermixed, from all participants).
3. `system-tests-compare.log`: The logs for the comparison to the reference results.

In addition, in the directories of the cases executed, you can find `system-tests-<case>.log` files.

### Numerical regressions

When the tests fail at the results comparison step, this typically means that there are numerical regressions (unless something wrong went undetected in a previous step). We use [fieldcompare](https://gitlab.com/dglaeser/fieldcompare) to compare all [preCICE exports](https://precice.org/configuration-export.html) to reference results generated from a previous run. Relevant files:

- `precice-exports/`: The coupling meshes of the test run.
- `reference-results/`: The coupling meshes of the reference run, as stored on Git LFS, expanded into `reference-results-unpacked`. For test cases using implicit coupling, the reference `.tar.gz` also contains the reference `precice-*-iterations.log` files.
- `diff-results/`: Numerical difference of the results in the two directories (computed with `fieldcompare dir --diff precice-exports/ reference/`). These are only present on failed comparisons.
- `iterations-logs/`: The `precice-*-iterations.log` files of the test run. Only present in test cases using implicit coupling. The comparisons to references only take into account the file SHA-256 checksums.

To reproduce the comparison locally, use the [same fieldcompare command](https://github.com/precice/tutorials/blob/develop/tools/tests/docker-compose.field_compare.template.yaml):

```bash
fieldcompare dir precice-exports/ reference-results-unpacked/<case>/ \
             --ignore-missing-reference-files \
             --ignore-unsupported-file-formats \
             -rtol 3e-7
```

The default relative tolerance (`-rtol`) is `3e-7`. Per-test overrides are possible in `tests.yaml` (e.g., `tolerance: 1e-2`). Set `skip_compare: true` to skip the comparison step and only verify that the build and run steps succeed.

The differences are only shown per file, and there is no global metric or other summary (see [related discussion in fieldcompare](https://gitlab.com/dglaeser/fieldcompare/-/work_items/69)).

Alternatively, [visualize the `precice-exports/diff_*.vtu` in ParaView](https://precice.org/configuration-export.html#visualization-with-paraview).

### Re-running from CI artifacts

When a system test fails in CI, download the **full** artifact:

`system_tests_run_<run_id>_<run_attempt>_full`

(a smaller `_logs` archive contains only log files). The archive contains a shared `runs/` directory:

```text
runs/
├── tools/                              # Dockerfiles and helpers (shared)
└── <tutorial>_<cases>_<timestamp>/     # one folder per system test
    ├── docker-compose.tutorial.yaml
    ├── docker-compose.field_compare.yaml   # written at build time when compare is configured
    ├── rerun-system-test.sh
    ├── system-tests-build.log
    ├── system-tests-run.log
    ├── system-tests-compare.log
    └── …
```

To re-run one test locally:

1. Extract the zip and keep the `runs/` layout (the test folder needs the sibling `tools/` directory).
2. `cd` into the test folder.
3. Run `./rerun-system-test.sh` (or `sh rerun-system-test.sh`).

The script rebuilds images, runs the tutorial, and (if present) runs fieldcompare with `--exit-code-from field-compare`, matching the Python runner. Compose paths are relative to the test folder (`..` is the parent `runs/` directory), so you can move the extracted tree elsewhere on a Linux host with Docker.

`docker-compose.field_compare.yaml` is written when the test is prepared for Docker build. The replay script fixes common permission issues from extracted archives.

Fieldcompare requires reference results in the artifact. If not already unpacked during the original CI run, unpack them manually first.

## Extending

### Adding new tests

Tests and test suites are defined in [`tests.yaml`](https://github.com/precice/tutorials/blob/develop/tools/tests/tests.yaml). By convention, every tutorial defines a test suite with the same name as its directory, and several test cases using combinations of the available participants. These test cases are later referenced by other test suites: these are typically the `release` and the test suites of different tested components.

The available cases are listed in the `metadata.yaml` of each tutorial. To add a new tutorial case as a test, add it to `metadata.yaml` and then define a test using it. Include that test in the relevant test suites.

Use the `max_time` or `max_time_windows` parameters to restrict the runtime of the test to the first few coupling time windows, to save time. Aim for a runtime of less than a minute (assuming cached components), if possible.

Some tutorials require setup before the simulation (e.g. switching configuration files). Use optional `run-before` and `run-after` fields in `tests.yaml` to run shell commands in the copied tutorial directory after copying and before Docker build (`run-before`), or after the simulation and before field comparison (`run-after`). Example:

```yaml
run-before: ./set-case.sh 1d3d
```

You will need to define a reference results file. The reference results can and should be generated on GitHub using the [Generate reference results (manual)](https://github.com/precice/tutorials/actions/workflows/generate-reference-results-manual.yml) workflow for the respective test suite. You might want to temporarily set the `selected` test suite for requesting results only for a subset of test cases.

When generating reference results for a feature branch, pass `TUTORIALS_REF` via `--build_args` (locally) or the optional `build_args` workflow input. The [Generate reference results (manual)](https://github.com/precice/tutorials/actions/workflows/generate-reference-results-manual.yml) workflow defaults to `TUTORIALS_REF:<from_ref>` when `build_args` is left empty, so you normally do not need to edit [`reference_versions.yaml`](reference_versions.yaml).

The results will be added to a Git LFS, but you will need special push access: just use the aforementioned GitHub Actions workflow, instead.

#### External test sources

Local test cases are defined in the `tutorials:` list. For test cases maintained in other git repositories or available as archives, use the `external:` list. A test suite may define both lists so that external cases can be referenced via YAML anchors alongside local tutorials (for example in the `release` suite).

```yaml
test_suites:
  some-external-test-case:
    external:
      - &some-external-test-case
        source:
          type: git
          url: https://github.com/some-user/tutorials.git
          ref: some-branch
          subdir: .
        path: path-to-test-case
        case_combination:
          ...

  release:
    tutorials:
      - *quickstart_openfoam_cpp
    external:
      - *some-external-test-case
```

Supported `source.type` values:

- `git`: shallow-clone `url` at `ref`. The `ref` must exist on the remote; clone failures are reported as errors.
- `archive`: download and extract a `.tar.gz` / `.zip` from `url`.

Fetched external sources are cached under `~/.cache/precice-tutorials` (or `PRECICE_EXTERNAL_CACHE_DIR`) so that repeated test runs and CI matrix jobs do not re-download the same repository or archive every time. The cache key is derived from the source URL, ref, and optional subdir.

The runner fetches the tutorial, copies it into the run directory, and then continues with the usual Docker build/run and fieldcompare steps. `TUTORIALS_REF` / `TUTORIALS_PR` build arguments apply only to test cases sourced from the tutorials repository; external test cases are pinned by `source.ref`. Reference result paths are resolved relative to the root directory of the fetched test case.

### Adding new components

To add a new component, a few changes are needed:

1. In the `components.yaml`, add a new component, defining all the parameters that it might need.
    - Add these parameters into the `reference_versions.yaml`.
    - Add fields for these parameters into the `system-tests-latest-components.yml` workflow.
2. In the `dockerfiles/<platform>/Dockerfile`, define a new stage for building your new component. Use the parameters you defined above.
    - Defining a `<component>_PR` variable will let you integrate the system tests into the respective component repository.
3. In the `component_templates.yaml`, define a component template. These are Jinja templates that are used by Docker Compose, and the main detail is how to run a simulation using that component.
4. Refer to the new component in the `metadata.yaml` of some tutorial and define some tests.

Some components might require injecting further changes to the configuration files of a tutorial, to make them use the testing setup (e.g., install specific versions of dependencies). You can set such steps for the `prepare:` service in the `docker-compose.template.yaml` file.

### Adding new repositories that can trigger the tests

If you want to trigger the system tests from a new repository:

1. Add a workflow file to that repository (under `.github/workflows/`). Example: [OpenFOAM adapter](https://github.com/precice/openfoam-adapter/blob/develop/.github/workflows/system-tests.yaml).
2. Create label as a trigger for workflow (under `Issues`->`Labels`->`New label`).
3. Give permissions: In the [preCICE organization settings](https://github.com/organizations/precice/settings), you need to add the new repository to the action secret `WORKFLOW_DISPATCH_TOKEN` and to the default actions runner group. Adding the new label directly to the pull request in which you add the workflow should already trigger the system tests.

## Implementation details

Each tutorial contains automation scripts (mainly `run.sh` and `clean.sh`), as well as metadata (`metadata.yaml`). The metadata file describes the available cases, how to run them, as well as their dependencies. A central `tests.yaml` file defines test suites, which execute different combinations of cases. The Python script `systemtests.py` executes the tests, filter for specific test suites.

Read more about the system tests in the publication [System Regression Tests for the preCICE Coupling Ecosystem](https://doi.org/10.14279/eceasst.v83.2614).

[![ECEASST](https://img.shields.io/badge/DOI-10.14279%2Feceasst.v83.2614-green)](https://doi.org/10.14279/eceasst.v83.2614)

<details markdown="1"><summary>Expand the implementation details...</summary>

### General architecture

Each tutorial directory contains a metadata file, describing which participants each case directory implements, and how to run it.

A list of tests describes all tests to be executed, grouped by test suites. Each test is a combination of tutorial cases.

Test steps include modifying the tutorial configuration files for the test system, building the Docker containers used by the respective Docker Compose service of each component, and comparing results to reference results using fieldcompare.

Tests are executed by the `systemtests.py` script, which starts the Docker Compose. This can be executed locally, and it is the same script that GitHub Actions also execute.

The multi-stage Docker build allows building each component separately from the same Dockerfile, while Docker reuses cached layers. The Docker Compose services consider GitHub Actions Cache when building the services, although the cache is currently only updated, but not hit (see https://github.com/precice/tutorials/pull/372#issuecomment-1748335750). For this reason, we are running on a dedicated self-hosted runner.

### File structure

Metadata and workflow/script files:

- `.github/workflows/`
  - `system-tests.yml`: workflow for running the tests, triggered by other workflows (e.g., other repositories)
  - `system-tests-latest-components.yml`: manual triggering front-end for `system-tests.yml`
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
  - `dockerfiles/`
    - Multi-stage build Dockerfiles that define how to build each component, in a layered approach
  - `docker-compose.template.yaml`: Describes how to prepare each test (Docker Compose service template)
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
  - `rerun-system-test.sh`: Helper script copied into each run directory so CI artifacts can be replayed locally.

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
- `path`: Where the tutorial is located, relative to the tutorials repository (or the tutorial folder name for external sources)
- `url`: A web page with more information on the tutorial
- `participants`: A list of preCICE participants, typically corresponding to different domains of the simulation
- `cases`: A list of solver configuration directories. Each element of the list includes:
  - `participant`: Which participant this solver case can serve as
  - `directory`: Where the case directory is located, relative to the tutorial directory
  - `run`: Command that executes the tutorial
  - `component`: Component or list of components that this case depends upon (typically an adapter)

### Components

The components mentioned in the Metadata are defined in the central `components.yaml` file. This file also specifies some arguments and their default values. These arguments can be anything, but most often they are version-related information. For example, the version of OpenFOAM used by the openfoam-adapter component, or the component version itself. For example:

```yaml
openfoam-adapter:
  template: component-templates/openfoam-adapter.yaml
  build_arguments:
    PLATFORM:
      default: "ubuntu_2404"
    PRECICE_REF:
      repository: https://github.com/precice/precice
      default: "develop"
    PRECICE_PRESET:
      default: "production-audit"
    TUTORIALS_REF:
      repository: https://github.com/precice/tutorials
      default: "develop"
    OPENFOAM_EXECUTABLE:
      default: "openfoam2512"
    OPENFOAM_ADAPTER_REF:
      repository: https://github.com/precice/openfoam-adapter
      default: "develop"
```

This `openfoam-adapter` component has the following attributes:

- `template`: A template for a Docker Compose service of this component
- `build_arguments`: Arguments passed to the Docker Compose service (arbitrary)

#### Naming schema for build_arguments

The following rules apply for the `build_arguments`:

- A build argument ending in `_REF` refers to a git commit-ish (like a tag or commit) being used to build the image.
- The `repository` parameter of any `_REF` argument points to the repository where the respective Git reference should be resolved. This is assumed to be consistent across components and the Dockerfile (and could be simplified).
- Some workflows set variables ending in `_PR`. These specify the GitHub pull request, which provides the above `_REF` and can be on a fork.
- A build argument ending in `_VERSION` refers to the version of a third-party dependency to use (e.g., DUNE).
- All other `build_arguments` are free of rules and up to the container maintainer.

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

- `image`: The base Docker image for this component, including a Git reference (tag), provided to the template as an argument (e.g., by the `systemtests.py` script).
- `depends_on`: Other services this service depends upon, typically a preparation service that fetches all components and tutorials.
- `volumes`: Directories mapped between the host and the container. Apart from directories relating to the users and groups, this also defines where to run the cases.
- `command`: How to run a case depending on this component, including how and where to redirect any screen output.

### Timeouts

The build and the run/compare steps use separate timeouts.

**Build:** Each component in `components.yaml` may set a `build_timeout` (default: 480s (8min)). The build step runs `docker compose build` once for all participant images, with one wall-clock subprocess timeout. That limit is the maximum `build_timeout` among the distinct components of the test (so that the slowest component is not cut off too early). You can override the default via the `PRECICE_SYSTEMTESTS_BUILD_TIMEOUT` environment variable.

**Run and compare:** Each test in `tests.yaml` may set a `timeout` (default: 180s (3min)), which applies to the running and results comparison steps only. You can override the default via the `PRECICE_SYSTEMTESTS_TIMEOUT` environment variable.

Tests can define a different `tolerance` in their `tests.yaml` entry, which applies to the fieldcompare step (relative tolerance, `-rtol`). The default is `3e-7`. Use `skip_compare: true` to skip fieldcompare entirely.

</details>
