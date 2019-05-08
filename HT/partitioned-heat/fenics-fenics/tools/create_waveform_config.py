from jinja2 import Environment, select_autoescape, FileSystemLoader
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-wr", "--waveform", nargs=2, default=[1, 1], type=int)
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

wr_tag = "-WR{N_Dirichlet}{N_Neumann}".format(N_Dirichlet=N_Dirichlet,
                                              N_Neumann=N_Neumann)

precice_config_name = "precice-config{wr_tag}.xml".format(wr_tag=wr_tag)

with open(precice_config_name, "w") as file:
    file.write(precice_config_template.render(temperatures=temperatures,
                                              fluxes=fluxes,
                                              convergence_limit=args.tolerance,
                                              total_time=total_time,
                                              window_size=window_size))

with open('precice-adapter-config-D{wr_tag}.json'.format(wr_tag=wr_tag), "w") as file:
    file.write(precice_adapter_D_template.render(N_Dirichlet=N_Dirichlet,
                                                 precice_config_name=precice_config_name))

with open('precice-adapter-config-N{wr_tag}.json'.format(wr_tag=wr_tag), "w") as file:
    file.write(precice_adapter_N_template.render(N_Neumann=N_Neumann,
                                                 precice_config_name=precice_config_name))
