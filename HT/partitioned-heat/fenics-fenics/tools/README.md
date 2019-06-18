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
