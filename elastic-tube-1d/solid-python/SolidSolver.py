from __future__ import division, print_function

import os
import sys
import argparse
import numpy as np
import precice
from precice import action_write_initial_data, action_read_iteration_checkpoint, \
    action_write_iteration_checkpoint


r0 = 1 / np.sqrt(np.pi)  # radius of the tube
a0 = r0**2 * np.pi  # cross sectional area
tau = 10**10  # timestep size, set it to a large value to enforce tau from precice_config.xml
N = 100  # number of elements in x direction
p0 = 0  # pressure at outlet
L = 10  # length of tube/simulation domain
E = 10000  # elasticity module
c_mk = np.sqrt(E / 2 / r0)  # wave speed


def crossSection0(N):
    return a0 * np.ones(N + 1)


print("Starting Solid Solver...")

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", nargs='?', type=str,
                    default="precice-config.xml")

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Did you forget adding the precice configuration file as an argument?")
    print("Try '$ python SolidSolver.py precice-config.xml'")
    quit()

print("N: " + str(N))

print("Configure preCICE...")
interface = precice.Interface("Solid", args.configurationFileName, 0, 1)
print("preCICE configured...")

dimensions = interface.get_dimensions()

pressure = p0 * np.ones(N + 1)
crossSectionLength = a0 * np.ones(N + 1)

meshID = interface.get_mesh_id("Solid-Nodes-Mesh")
crossSectionLengthID = interface.get_data_id("CrossSectionLength", meshID)
pressureID = interface.get_data_id("Pressure", meshID)

vertexIDs = np.zeros(N + 1)
grid = np.zeros([N + 1, dimensions])

grid[:, 0] = np.linspace(0, L, N + 1)  # x component
grid[:, 1] = 0  # np.linspace(0, config.L, N+1)  # y component, leave blank

vertexIDs = interface.set_mesh_vertices(meshID, grid)

t = 0

print("Solid: init precice...")

# preCICE defines timestep size of solver via precice-config.xml
precice_dt = interface.initialize()

if interface.is_action_required(action_write_initial_data()):
    interface.write_block_scalar_data(crossSectionLengthID, vertexIDs, crossSectionLength)
    interface.mark_action_fulfilled(action_write_initial_data())

interface.initialize_data()

if interface.is_read_data_available():
    pressure = interface.read_block_scalar_data(pressureID, vertexIDs)

crossSection0 = crossSection0(pressure.shape[0] - 1)
pressure0 = p0 * np.ones_like(pressure)

while interface.is_coupling_ongoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.is_action_required(action_write_iteration_checkpoint()):
        interface.mark_action_fulfilled(action_write_iteration_checkpoint())

    crossSectionLength = crossSection0 * (
        (pressure0 - 2.0 * c_mk ** 2) ** 2 / (pressure - 2.0 * c_mk ** 2) ** 2)

    interface.write_block_scalar_data(crossSectionLengthID, vertexIDs, crossSectionLength)
    precice_dt = interface.advance(precice_dt)
    pressure = interface.read_block_scalar_data(pressureID, vertexIDs)

    if interface.is_action_required(action_read_iteration_checkpoint()):  # i.e. not yet converged
        interface.mark_action_fulfilled(action_read_iteration_checkpoint())
    else:
        t += precice_dt

print("Exiting SolidSolver")

interface.finalize()
