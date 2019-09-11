from jinja2 import Environment, FileSystemLoader
import os, stat
from coupling_schemes import CouplingScheme
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=0.0, type=float)
parser.add_argument("-stol", "--solver-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**-12, type=float)
parser.add_argument("-dd", "--domain-decomposition", help="set kind of domain decomposition being used", default="DN", type=str, choices=['DN', 'ND'])
parser.add_argument("-subs", "--plain-subcycling", help="if set, do not interpolate between samples, but use plain subcycling for coupling", dest='plain_subcycling', action='store_true')

args = parser.parse_args()

if args.plain_subcycling:
    use_subcycling = "--plain-subcycling"
else:
    use_subcycling = ""

wr_lefts = [1, 2, 5]
wr_rights = [1, 2, 5]
window_sizes = [1.0, 0.5, 0.2, 0.1]
#coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

env = Environment(
    loader=FileSystemLoader('./templates')
)

configs_template = env.get_template('config_creation.sh')

with open('config_creation.sh', "w") as file:
    file.write(configs_template.render(window_sizes=window_sizes,
                                       wr_lefts=wr_lefts,
                                       wr_rights=wr_rights,
                                       coupling_schemes=coupling_schemes,
				       gamma=args.gamma,
                                       solver_tolerance=args.solver_tolerance,
                                       domain_decomposition=args.domain_decomposition,
                                       plain_subcycling=use_subcycling))

st = os.stat('config_creation.sh')
os.chmod('config_creation.sh', st.st_mode | stat.S_IEXEC)
