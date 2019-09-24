from jinja2 import Environment, FileSystemLoader
import os, stat
from coupling_schemes import CouplingScheme

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-wrl", "--wr-lefts", nargs='+', help="Choose WR setups on left half of domain", default=[1, 2, 3, 5])
parser.add_argument("-wrr", "--wr-rights", nargs='+', help="Choose WR setups on left half of domain", default=[1, 2, 3, 5])
parser.add_argument("-Dts", "--window-sizes", nargs='+', help="Choose Window Sizes being computed", default=[5.0, 2.0, 1.0, 0.5, 0.2])
args = parser.parse_args()

wr_lefts = args.wr_lefts
wr_rights = args.wr_rights
window_sizes = args.window_sizes
#coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

env = Environment(
    loader=FileSystemLoader('./templates')
)

configs_template = env.get_template('runexperiments.sh')

with open('runexperiments.sh', "w") as file:
    file.write(configs_template.render(window_sizes=window_sizes,
                                       wr_lefts=wr_lefts,
                                       wr_rights=wr_rights,
                                       coupling_schemes=coupling_schemes))

st = os.stat('runexperiments.sh')
os.chmod('runexperiments.sh', st.st_mode | stat.S_IEXEC)
