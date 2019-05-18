from jinja2 import Environment, select_autoescape, FileSystemLoader
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument("-wr", "--waveform", nargs=2, default=[1, 1], type=int)
parser.add_argument("-dT", "--window-size", default=1, type=float)
parser.add_argument("-tol", "--tolerance", default='1e-12', type=str)
args = parser.parse_args()

temperatures = []
fluxes = []

N_Dirichlet = args.waveform[0]
N_Neumann = args.waveform[1]

"""
define timestepping setup. Be aware of the following relationships:

1. N_Neumann * dt_Neumann = N_Dirichlet * dt_Dirichlet = window_size
2. window_size * N_coupling = total_time
"""
total_time = 1
dt = .1
window_size = dt * np.max([N_Neumann, N_Dirichlet])

for i in range(N_Neumann):
    temperatures.append("Temperature{i}".format(i=i+1))

for i in range(N_Dirichlet):
    fluxes.append("Flux{i}".format(i=i+1))

env = Environment(
    loader=FileSystemLoader('./templates'),
    autoescape=select_autoescape(['xml', 'json'])
)
precice_config_template = env.get_template('precice-config.xml')
precice_adapter_D_template = env.get_template('precice-adapter-config-D.json')
precice_adapter_N_template = env.get_template('precice-adapter-config-N.json')

wr_tag = "WR{N_Dirichlet}{N_Neumann}".format(N_Dirichlet=N_Dirichlet,
                                              N_Neumann=N_Neumann)
window_size = "dT{dT}".format(dT=args.window_size)

target_path = os.path.join("experiments", wr_tag, window_size)

if not os.path.exists(target_path):
    os.makedirs(target_path)

precice_config_name = "precice-config.xml"

with open(os.path.join( target_path, precice_config_name), "w") as file:
    file.write(precice_config_template.render(temperatures=temperatures,
                                              fluxes=fluxes,
                                              convergence_limit=args.tolerance,
                                              total_time=total_time,
                                              window_size=window_size))

with open(os.path.join( target_path, 'precice-adapter-config-D.json'), "w") as file:
    file.write(precice_adapter_D_template.render(N_Dirichlet=N_Dirichlet,
                                                 precice_config_name=precice_config_name))

with open(os.path.join( target_path, 'precice-adapter-config-N.json'), "w") as file:
    file.write(precice_adapter_N_template.render(N_Neumann=N_Neumann,
                                                 precice_config_name=precice_config_name))
