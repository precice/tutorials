from jinja2 import Environment, select_autoescape, FileSystemLoader
import argparse
import numpy as np
import os, stat
from participants import Participant

wr_lefts = [1,2,3,5,10]
wr_rights = [1,2,3,5,10]
window_sizes = [1.0, 0.5, 0.2, 0.1]
first_participants = [Participant.DIRICHLET.name, Participant.NEUMANN.name]

env = Environment(
    loader=FileSystemLoader('./templates')
)

configs_template = env.get_template('config_creation.sh')

with open('config_creation.sh', "w") as file:
    file.write(configs_template.render(window_sizes=window_sizes,
                                       wr_lefts=wr_lefts,
                                       wr_rights=wr_rights,
                                       first_participants=first_participants))

st = os.stat('config_creation.sh')
os.chmod('config_creation.sh', st.st_mode | stat.S_IEXEC)
