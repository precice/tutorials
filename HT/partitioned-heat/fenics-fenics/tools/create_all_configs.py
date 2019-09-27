from jinja2 import Environment, FileSystemLoader
import os, stat
from coupling_schemes import CouplingScheme
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-stol", "--solver-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**-12, type=float)
parser.add_argument("-qntol", "--quasi-newton-tolerance", help="set accepted error in the quasi newton scheme", default='1e-12', type=str)
parser.add_argument("-dd", "--domain-decomposition", help="set kind of domain decomposition being used", default="DN", type=str, choices=['DN', 'ND'])
parser.add_argument("-t", "--time-dependence", help="choose whether there is a linear (l), quadratic (q) or sinusoidal (s) dependence on time", type=str, default="l")
parser.add_argument("-mth", "--method", help="time stepping method being used", default='ie')
parser.add_argument("-wri", "--waveform-interpolation-strategy", help="specify interpolation strategy used by waveform relaxation", default="linear", choices=['linear', 'quadratic', 'cubic', 'quartic'], type=str)
parser.add_argument("-exec", "--executable", help="choose name of executable", default='heat.py')
parser.add_argument("-pp", "--post-processing", help="specify postprocessing scheme used by preCICE", default="qn-active", choices=['none', 'underrelaxation', 'qn-passive', 'qn-passive-fair', 'qn-active', 'qn-active-fair'], type=str)
parser.add_argument("-m", "--monolithic", dest='monolithic', help="switch to monolithic case", action='store_true')
parser.add_argument("-p", "--partitioned", dest='monolithic', help="switch to partitioned case", action='store_false')
parser.add_argument("-T", "--simulation-time", default=10, type=float)
parser.add_argument("-Dts", "--window-sizes", nargs='+', help="Choose Window Sizes being computed", default=[5.0, 2.0, 1.0, 0.5, 0.2, 0.1])
parser.add_argument("-wrl", "--wr-lefts", nargs='+', help="Choose WR setups on left half of domain", default=[1, 2, 3, 5])
parser.add_argument("-wrr", "--wr-rights", nargs='+', help="Choose WR setups on left half of domain", default=[1, 2, 3, 5])

args = parser.parse_args()

wr_lefts = args.wr_lefts
wr_rights = args.wr_rights
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

env = Environment(
    loader=FileSystemLoader('./templates')
)

configs_template = env.get_template('config_creation.sh')

if args.monolithic:
    case_flag = '--monolithic'
else:
    case_flag = '--partitioned'

configuration = vars(args)

import json

with open('configuration.json', 'w') as fp:
    json.dump(configuration, fp)

with open('config_creation.sh', "w") as file:
    file.write(configs_template.render(window_sizes=args.window_sizes,
                                       wr_lefts=wr_lefts,
                                       wr_rights=wr_rights,
                                       coupling_schemes=coupling_schemes,
                                       solver_tolerance=args.solver_tolerance,
                                       domain_decomposition=args.domain_decomposition,
                                       time_dependence=args.time_dependence,
                                       method=args.method,
                                       executable=args.executable,
                                       waveform_interpolation_strategy=args.waveform_interpolation_strategy,
                                       quasi_newton_tolerance=args.quasi_newton_tolerance,
                                       post_processing=args.post_processing,
                                       case_flag=case_flag,
                                       simulation_time=args.simulation_time))

st = os.stat('config_creation.sh')
os.chmod('config_creation.sh', st.st_mode | stat.S_IEXEC)
