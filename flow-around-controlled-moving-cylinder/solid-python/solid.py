from __future__ import division

import argparse
import numpy as np
import precice

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", default="../precice-config.xml",
                    help="Name of the xml config file.", type=str)

try:
    args = parser.parse_args()
    print(args)
except SystemExit:
    print("")
    print("Usage: python3 ./solverdummy.py precice-config.xml")
    quit()

# mass-spring-damper system parameters
mass = 0.03575  # kg
k_spring = 69.48  # N/m
d_damper = 0.0043  # N s/m

state_old = np.zeros(3)
state = np.zeros(3)


def update(m, k, d, state, dt, force, controlled_spring_displacement):
    x_old, v_old, a_old = state
    x = x_old + v_old * dt + 0.5 * a_old * dt**2
    a = (force - d * v_old - k * x + k * controlled_spring_displacement) / m
    v = v_old + 0.5 * (a_old + a) * dt
    return np.array([x, v, a])


configuration_file_name = args.configurationFileName

participant_name = "Solid"
write_data_name = "Displacement-Cylinder"
read_data_name_force = "Force"
read_data_name_displacement = "Displacement-Spring"
mesh_name = "Mesh-Solid"

num_vertices = 1  # Number of vertices
solver_process_index = 0
solver_process_size = 1

participant = precice.Participant(participant_name, configuration_file_name,
                                  solver_process_index, solver_process_size)

assert (not participant.requires_mesh_connectivity_for(mesh_name))

dimensions = participant.get_mesh_dimensions(mesh_name)

vertices = np.zeros((num_vertices, dimensions))
read_data_force = np.zeros((num_vertices, dimensions))
read_data_displacement = np.zeros((num_vertices, dimensions))
write_data = np.zeros((num_vertices, dimensions))

vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

# initialize data
if participant.requires_initial_data():
    participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

participant.initialize()
t = 0

while participant.is_coupling_ongoing():
    if participant.requires_writing_checkpoint():
        print("Writing checkpoint")
        state_old_cp = state_old
        state_cp = state

    dt = participant.get_max_time_step_size()

    read_data_force = participant.read_data(mesh_name, read_data_name_force, vertex_ids, dt)
    force = read_data_force[0, 1]
    read_data_displacement = participant.read_data(mesh_name, read_data_name_displacement, vertex_ids, dt)
    displacement_spring = read_data_displacement[0, 1]

    # compute next time step
    state = update(mass, k_spring, d_damper, state_old, dt, force, displacement_spring)

    # cylinder is fixed in x-direction
    write_data[0, 0] = 0

    # cylinder moves in y-direction according to spring-damper-mass-equation
    write_data[0, 1] = state[0]

    participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

    print("Solid: Advancing in time")
    participant.advance(dt)

    if participant.requires_reading_checkpoint():
        print("Reading checkpoint")
        state_old = state_old_cp
        state = state_cp

    else:
        state_old = state
        t = t + dt

participant.finalize()
