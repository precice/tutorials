from jinja2 import Environment, FileSystemLoader
import os, stat
from coupling_schemes import CouplingScheme
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=0.0, type=float)
parser.add_argument("-stol", "--solver-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**-12, type=float)

args = parser.parse_args()                                    

wr_lefts = [1, 2, 5, 10]
wr_rights = [1, 2, 5, 10]
window_sizes = [1.0, 0.5, 0.2, 0.1]
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]

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
                                       solver_tolerance=args.solver_tolerance))

st = os.stat('config_creation.sh')
os.chmod('config_creation.sh', st.st_mode | stat.S_IEXEC)
