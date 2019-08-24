from __future__ import division

import argparse
import numpy as np
import precice_future as precice

force_excitement = np.array([0, 10, 0])
omega_excitement = 1

# beam geometry
L = 0.35  # length
H = 0.02  # height
y_bottom = 0.2 - 0.5 * H  # y coordinate of bottom surface of beam
y_top = y_bottom + H  # y coordinate of top surface of beam
x_left = 0.25  # x coordinate of left surface of beam
x_right = x_left + L  # x coordinate of right surface of beam
x_attack = x_right
y_attack = .5 * (y_bottom + y_top)


def compute_force(time):
    omega = omega_excitement
    force_0 = force_excitement
    return force_0 * np.sin(omega * time)


parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy precice-config participant-name mesh-name")
    quit()

configuration_file_name = args.configurationFileName
participant_name = "Excitement"
mesh_name = "Excitement-Mesh"
write_data_name = "Forces1"
read_data_name = "Displacements1"

n = 1

solver_process_index = 0
solver_process_size = 1

interface = precice.Interface(participant_name, solver_process_index, solver_process_size)
interface.configure(configuration_file_name)

mesh_id = interface.get_mesh_id(mesh_name)

dimensions = interface.get_dimensions()
vertices = np.array([x_attack, y_attack, 0])

vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)
write_data_id = interface.get_data_id(write_data_name, mesh_id)
read_data_id = interface.get_data_id(read_data_name, mesh_id)

dt = interface.initialize()
t = 0

if interface.is_action_required(precice.action_write_initial_data()):
    forces = compute_force(t + dt)
    print("write initial forces: {}".format(forces))
    interface.write_block_vector_data(write_data_id, vertex_ids, forces)
    interface.fulfilled_action(precice.action_write_initial_data())

interface.initialize_data()

while interface.is_coupling_ongoing():

    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        print("DUMMY: Writing iteration checkpoint")
        interface.fulfilled_action(precice.action_write_iteration_checkpoint())

    forces = compute_force(t+dt)
    print("sending forces: {}".format(forces))
    interface.write_block_vector_data(write_data_id, vertex_ids, forces)
    dt = interface.advance(dt)
    displacement = interface.read_block_vector_data(read_data_id, vertex_ids)

    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        print("DUMMY: Reading iteration checkpoint")
        interface.fulfilled_action(precice.action_read_iteration_checkpoint())
    else:
        t += dt
        print("DUMMY: Advancing in time")

interface.finalize()
print("DUMMY: Closing python solver dummy...")

