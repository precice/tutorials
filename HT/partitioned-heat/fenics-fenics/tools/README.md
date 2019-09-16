# Configuration Generator

The enclosed python script allows you to generate configuration files for different waveform relaxation setups.

## Dependencies

* `jinja2`: `pip3 install --user jinja2`

## Running

Run `python3 create_waveform_config.py -wr 1 2` for creating the configuration files that you need for the `WR12` setup. This includes:

* `precice-config.xml`
* `precice-adapter-config-D.json`
* `precice-adapter-config-N.json`

The templates that are used for generating these files can be found in `templates`. The generated files are created in `experiments/WRXX/dT00/*`.

## Create the whole cross product of configs

### Set of parameters

There are many different parameters for our experiments:

* Number of substeps left, `wr_lefts = [1,2,3,5,10]`
* Number of substeps right, `wr_rights = [1,2,3,5,10]`
* Size of the window, `window_sizes = [1.0, 0.5, 0.2, 0.1]`
* Coupling schemes, `coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]`

Other parameters are fixed. E.g.:

* total simulation time `T=10`

### How to create all configs

*Note:* `$ROOT` denotes `tutorials/partitionced-heat/fenics-fenics/`

1. running `python3 create_all_configs.py` will create a file `config_creation.sh`. Optional parameters for setting `gamma` and/or `solver_tolerance` may be used.
2. running `./config_creation.sh` will execute `python3 create_waveform_config.py` for all possible configurations and will create all configs in the `experiments` folder
3. copy the `experiments` folder to `tutorials/partitionced-heat/fenics-fenics/`
4. running `python3 create_runexperiments.py` will create a file `runexperiments.sh`
5. copy `runexperiments.sh` to `$ROOT`
6. running `$ROOT/runexperiments.sh` will run all experiments with the configuration files in `$ROOT/experiments`

### Postprocessing

Run `python3 postproc.py` from **this folder**

### Create a report

Go to parent folder of this one and run the following commands:

```
cd ..
python3 tools/create_report.py -p experiments
pandoc report.md -o report.pdf
```

*Remark* use the `-p` flag to provide the path to the folder containing the experiments.

