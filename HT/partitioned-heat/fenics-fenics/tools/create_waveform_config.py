from jinja2 import Environment, select_autoescape, FileSystemLoader
import argparse
import os, stat
from coupling_schemes import CouplingScheme


parser = argparse.ArgumentParser()
parser.add_argument("-wr", "--waveform", nargs=2, default=[1, 1], type=int)
parser.add_argument("-dT", "--window-size", default=1.0, type=float)
parser.add_argument("-T", "--simulation-time", default=10, type=float)
parser.add_argument("-qntol", "--quasi-newton-tolerance", help="set accepted error in the quasi newton scheme", default='1e-12', type=str)
parser.add_argument("-cpl", "--coupling-scheme", default=CouplingScheme.SERIAL_FIRST_DIRICHLET.name, type=str)
parser.add_argument("-g", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=0.0, type=float)
parser.add_argument("-stol", "--solver-tolerance", help="set accepted error of numerical solution w.r.t analytical solution", default=10**-12, type=float)
args = parser.parse_args()

temperatures = []
fluxes = []

N_Dirichlet = args.waveform[0]
N_Neumann = args.waveform[1]
if args.coupling_scheme == CouplingScheme.SERIAL_FIRST_DIRICHLET.name:
    coupling_scheme = CouplingScheme.SERIAL_FIRST_DIRICHLET
elif args.coupling_scheme == CouplingScheme.SERIAL_FIRST_NEUMANN.name:
    coupling_scheme = CouplingScheme.SERIAL_FIRST_NEUMANN
elif args.coupling_scheme == CouplingScheme.PARALLEL.name:
    coupling_scheme = CouplingScheme.PARALLEL
else:
    raise Exception("invalid input {} for --coupling-scheme".format(args.coupling_scheme))

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

if coupling_scheme == CouplingScheme.SERIAL_FIRST_DIRICHLET:
    precice_config_template = env.get_template('precice-config_serialImplicit_firstDirichlet.xml')
elif coupling_scheme == CouplingScheme.SERIAL_FIRST_NEUMANN:
    precice_config_template = env.get_template('precice-config_serialImplicit_firstNeumann.xml')
elif coupling_scheme == CouplingScheme.PARALLEL:
    precice_config_template = env.get_template('precice-config_parallelImplicit.xml')

precice_adapter_D_template = env.get_template('precice-adapter-config-D.json')
precice_adapter_N_template = env.get_template('precice-adapter-config-N.json')
runall_template = env.get_template("runall.sh")

wr_tag = "WR{N_Dirichlet}{N_Neumann}".format(N_Dirichlet=N_Dirichlet,
                                             N_Neumann=N_Neumann)
window_tag = "dT{dT}".format(dT=args.window_size)
coupling_tag = "{}".format(coupling_scheme.name)
total_time = args.simulation_time
target_path = os.path.join("experiments", wr_tag, window_tag, coupling_tag)

if not os.path.exists(target_path):
    os.makedirs(target_path)

precice_config_name = "precice-config.xml"

with open(os.path.join( target_path, precice_config_name), "w") as file:
    file.write(precice_config_template.render(temperatures=temperatures,
                                              fluxes=fluxes,
                                              convergence_limit=args.quasi_newton_tolerance,
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
    file.write(runall_template.render(wr_dirichlet=N_Dirichlet,
                                      wr_neumann=N_Neumann,
                                      window_size=args.window_size,
                                      coupling_scheme=coupling_scheme.name,
				      gamma=args.gamma,
                                      error_tolerance=args.solver_tolerance))

st = os.stat(runall_path)
os.chmod(runall_path, st.st_mode | stat.S_IEXEC)
