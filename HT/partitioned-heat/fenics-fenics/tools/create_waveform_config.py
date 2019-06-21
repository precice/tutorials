from jinja2 import Environment, select_autoescape, FileSystemLoader
import argparse
import numpy as np
import os, stat
from enum import Enum
from participants import Participant


parser = argparse.ArgumentParser()
parser.add_argument("-wr", "--waveform", nargs=2, default=[1, 1], type=int)
parser.add_argument("-dT", "--window-size", default=1.0, type=float)
parser.add_argument("-T", "--simulation-time", default=10, type=float)
parser.add_argument("-tol", "--tolerance", default='1e-12', type=str)
parser.add_argument("-first", "--first-participant", default="DIRICHLET", type=str)
args = parser.parse_args()

temperatures = []
fluxes = []

N_Dirichlet = args.waveform[0]
N_Neumann = args.waveform[1]
if args.first_participant == Participant.DIRICHLET.name:
    first_participant = Participant.DIRICHLET
elif args.first_participant == Participant.NEUMANN.name:
    first_participant = Participant.NEUMANN
else:
    raise Exception("invalid input {} for --first-participant".format(args.first_participant))

"""
define timestepping setup. Be aware of the following relationships:

1. N_Neumann * dt_Neumann = N_Dirichlet * dt_Dirichlet = window_size
2. window_size * N_coupling = total_time
"""

for i in range(N_Neumann):
    temperatures.append("Temperature{i}".format(i=i+1))

for i in range(N_Dirichlet):
    fluxes.append("Flux{i}".format(i=i+1))

env = Environment(
    loader=FileSystemLoader('./templates'),
    autoescape=select_autoescape(['xml', 'json'])
)

if first_participant == Participant.DIRICHLET:
    precice_config_template = env.get_template('precice-config_firstDirichlet.xml')
elif first_participant == Participant.NEUMANN:
    precice_config_template = env.get_template('precice-config_firstNeumann.xml')

precice_adapter_D_template = env.get_template('precice-adapter-config-D.json')
precice_adapter_N_template = env.get_template('precice-adapter-config-N.json')
runall_template = env.get_template("runall.sh")

wr_tag = "WR{N_Dirichlet}{N_Neumann}".format(N_Dirichlet=N_Dirichlet,
                                             N_Neumann=N_Neumann)
window_tag = "dT{dT}".format(dT=args.window_size)
first_participant_tag = "first_{}".format(first_participant.name)
total_time = args.simulation_time
target_path = os.path.join("experiments", wr_tag, window_tag, first_participant_tag)

if not os.path.exists(target_path):
    os.makedirs(target_path)

precice_config_name = "precice-config.xml"

with open(os.path.join( target_path, precice_config_name), "w") as file:
    file.write(precice_config_template.render(temperatures=temperatures,
                                              fluxes=fluxes,
                                              convergence_limit=args.tolerance,
                                              total_time=total_time,
                                              window_size=args.window_size))

with open(os.path.join( target_path, 'precice-adapter-config-D.json'), "w") as file:
    file.write(precice_adapter_D_template.render(N_Dirichlet=N_Dirichlet,
                                                 precice_config_name=precice_config_name))

with open(os.path.join( target_path, 'precice-adapter-config-N.json'), "w") as file:
    file.write(precice_adapter_N_template.render(N_Neumann=N_Neumann,
                                                 precice_config_name=precice_config_name))

runall_path = os.path.join( target_path, 'runall.sh')
with open(runall_path, "w") as file:
    file.write(runall_template.render(wr_left=N_Dirichlet,
                                      wr_right=N_Neumann,
                                      window_size=args.window_size,
                                      first_participant=first_participant.name))

st = os.stat(runall_path)
os.chmod(runall_path, st.st_mode | stat.S_IEXEC)
