from __future__ import division

import argparse
import numpy as np
import precice

configuration_file_name = "../precice-config.xml"
participant_name = "Fluid"
mesh_name = "Fluid-Mesh"
write_data_name = 'Force'

solver_process_index = 0
solver_process_size = 1

interface = precice.Participant(participant_name, configuration_file_name, solver_process_index, solver_process_size)
dimensions = interface.get_mesh_dimensions(mesh_name)
assert (dimensions == 2)

# infomation about beam geometry and position
W = 0.1  # width of the beam
H = 1.0  # height of the beam
x_left = 0.0 - 0.5 * W  # x-coordinate of the left side of the beam
y_bottom = 0.0  # bottom of the beam
y_top = y_bottom + H  # top of the beam
F_max = 1.0
# define vertices where we want to prove a load
n = 100  # Number of vertices
vertices = np.zeros((n, dimensions))
vertices[:, 0] = x_left  # all vertices are at left side of beam
vertices[:, 1] = np.linspace(y_bottom, y_top, n)  # have n equally disrtibuted vertices
vertex_ids = interface.set_mesh_vertices(mesh_name, vertices)

interface.initialize()
solver_dt = np.inf  # we just want to use dt = precice_dt

while interface.is_coupling_ongoing():
    if interface.requires_writing_checkpoint():
        pass

    precice_dt = interface.get_max_time_step_size()
    dt = min([solver_dt, precice_dt])

    write_data = np.zeros((n, dimensions))
    write_data[:, 0] = F_max * vertices[:, 1] / H  # linearly increasing load
    write_data[:, 1] = 0

    interface.write_data(mesh_name, write_data_name, vertex_ids, write_data)
    interface.advance(dt)

    if interface.requires_reading_checkpoint():
        pass

interface.finalize()
