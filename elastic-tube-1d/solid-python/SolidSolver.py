from __future__ import division, print_function

import os
import sys
import argparse
import numpy as np
import precice

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
interface = precice.Participant("Solid", args.configurationFileName, 0, 1)
print("preCICE configured...")

meshName = "Solid-Nodes-Mesh"
crossSectionLengthName = "CrossSectionLength"
pressureName = "Pressure"

dimensions = interface.get_mesh_dimensions(meshName)

pressure = p0 * np.ones(N + 1)
crossSectionLength = a0 * np.ones(N + 1)

# Define mesh coordinates and register coordinates
grid = np.zeros([N + 1, dimensions])
grid[:, 0] = np.linspace(0, L, N + 1)  # x component
grid[:, 1] = 0  # np.linspace(0, config.L, N+1)  # y component, leave blank

vertexIDs = interface.set_mesh_vertices(meshName, grid)

if interface.requires_initial_data():
    interface.write_data(meshName, crossSectionLengthName,
                         vertexIDs, crossSectionLength)

print("Solid: init precice...")

interface.initialize()

# Calculate initial crossSection and local pressure from initial received pressure
pressure = interface.read_data(meshName, pressureName, vertexIDs, 0)

crossSection0 = crossSection0(pressure.shape[0] - 1)
pressure0 = p0 * np.ones_like(pressure)


def solve(pressure):
    return crossSection0 * ((pressure0 - 2.0 * c_mk ** 2) ** 2 / (pressure - 2.0 * c_mk ** 2) ** 2)


while interface.is_coupling_ongoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.requires_writing_checkpoint():
        pass

    precice_dt = interface.get_max_time_step_size()

    # Read data, solve timestep, and write data
    pressure = interface.read_data(
        meshName, pressureName, vertexIDs, precice_dt)
    crossSectionLength = solve(pressure)
    interface.write_data(meshName, crossSectionLengthName,
                         vertexIDs, crossSectionLength)

    # Advance the coupling
    interface.advance(precice_dt)

    if interface.requires_reading_checkpoint():
        pass

print("Exiting SolidSolver")
