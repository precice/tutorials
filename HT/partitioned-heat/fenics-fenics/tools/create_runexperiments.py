from jinja2 import Environment, FileSystemLoader
import os, stat
from coupling_schemes import CouplingScheme

wr_lefts = [1, 2, 5]
wr_rights = [1, 2, 5]
window_sizes = [1.0]#, 0.5, 0.2, 0.1]
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
