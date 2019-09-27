from jinja2 import Environment, FileSystemLoader
import os, stat
from tools.coupling_schemes import CouplingScheme
from problem_setup import get_manufactured_solution
from tools.postproc_errors import create_error_table
from tools.postproc import create_qn_table
import argparse
import sympy as sp
import datetime
import codecs
import json

parser = argparse.ArgumentParser()
parser.add_argument('--prefix', '-p', help='Path prefix for results', type=str, default=os.path.join(os.path.dirname(os.path.realpath(__file__)),'experiments'))
parser.add_argument('--code-prefix', '-cp', help='Path prefix for code', type=str, default='~')
parser.add_argument("--configuration", type=str, default="experiments/configuration.json")

args = parser.parse_args()

with open(args.configuration, 'r') as fp:
    config = json.load(fp)

print(config)

alpha, beta, gamma = 3.0, 1.2, 1.0

manufactured_solution = sp.latex(get_manufactured_solution(config['time_dependence'], alpha, beta, gamma).subs("x[0]","x").subs("x[1]","y"))

env = Environment(
    loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)),'tools', 'templates'))
)

experiment_path = os.path.realpath(args.prefix)

if config['monolithic']:
    case_flag = '--monolithic'
else:
    case_flag = '--partitioned'

evaluated_wr = []
for wr_l in config['wr_lefts']:
    for wr_r in config['wr_rights']:
        evaluated_wr.append(int("{}{}".format(wr_l, wr_r)))
evaluated_dT = config['window_sizes']
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

print("parse experiments from {}".format(experiment_path))

error_table, _, _ = create_error_table(experiment_path, evaluated_wr, evaluated_dT, coupling_schemes, config['simulation_time'])

qn_table, _, _ = create_qn_table(experiment_path, evaluated_wr, evaluated_dT, coupling_schemes, config['simulation_time'])

report_template = env.get_template('report.md')

import subprocess
adapter_path = os.path.expanduser(os.path.join(args.code_prefix, "fenics-adapter"))
precice_path = os.path.expanduser(os.path.join(args.code_prefix, "precice"))
waveform_bindings_path = os.path.expanduser(os.path.join(args.code_prefix, "waveform-bindings"))

tutorials_hash = subprocess.check_output(["git", "describe", "--always", "--dirty"]).strip().decode() + " from " + subprocess.check_output(["git", "ls-remote", "--get-url"]).strip().decode()
adapter_hash = subprocess.check_output(["git", "-C", adapter_path, "describe", "--always", "--dirty"]).strip().decode() + " from " + subprocess.check_output(["git", "-C", adapter_path, "ls-remote", "--get-url"]).strip().decode()
if config['executable'] == "heat.py":
    waveform_bindings_hash = subprocess.check_output(["git", "-C", waveform_bindings_path, "describe", "--always", "--dirty"]).strip().decode() + " from " + subprocess.check_output(["git", "-C", waveform_bindings_path, "ls-remote", "--get-url"]).strip().decode()
else:
    waveform_bindings_hash = "/"
precice_hash = subprocess.check_output(["git", "-C", precice_path, "describe", "--always", "--dirty"]).strip().decode() + " from " + subprocess.check_output(["git", "-C", precice_path, "ls-remote", "--get-url"]).strip().decode()

try:
    with codecs.open('report.md', "w", 'utf-8') as file:
        file.write(report_template.render(date=str(datetime.datetime.now()),
                                      alpha=alpha,
                                      beta=beta,
                                      gamma=gamma,
                                      manufactured_solution=manufactured_solution,
                                      qn_table=qn_table,
                                      error_table=error_table,
                                      tutorials_hash=tutorials_hash,
                                      adapter_hash=adapter_hash,
                                      waveform_bindings_hash=waveform_bindings_hash,
                                      precice_hash=precice_hash,
                                      method=config['method'],
                                      executable=config['executable'].replace('_','\_'),
                                      time_dependence=config['time_dependence'],
                                      waveform_interpolation_strategy=config['waveform_interpolation_strategy'],
                                      coupling_tolerance=config['quasi_newton_tolerance'],
                                      solver_tolerance=config['solver_tolerance'],
                                      post_processing=config['post_processing'],
                                      simulation_time=config['simulation_time'],
                                      case_flag=case_flag))

except UnicodeEncodeError:
    raise UnicodeEncodeError
