from __future__ import division

import argparse
import numpy as np
import precice
from matplotlib import pyplot as plt

configuration_file_name = "../precice-config.xml"
participant_name = "Fluid"
mesh_name = "Fluid-Mesh"
write_data_name = 'Force'

num_vertices = 100  # Number of vertices

solver_process_index = 0
solver_process_size = 1

interface = precice.Participant(participant_name, configuration_file_name, solver_process_index, solver_process_size)

dimensions = interface.get_mesh_dimensions(mesh_name)
vertices = np.zeros((num_vertices, dimensions))
write_data = np.zeros((num_vertices, dimensions))

W = .1
H = 1
x_left = 0 - W / 2
y_bottom = 0
y_top = y_bottom + H
F_max = 1

vertices[:, 0] = x_left  # all vertices are at left side of beam
vertices[:, 1] = np.linspace(y_bottom, y_top, num_vertices)  # have num_vertices equally disrtibuted vertices

write_data[:, 0] = F_max * vertices[:, 1] / H  # linearly increasing load
write_data[:, 1] = 0

vertex_ids = interface.set_mesh_vertices(mesh_name, vertices)

interface.initialize()

t = 0
solver_dt = np.inf  # we just want to use dt = precice_dt
time = []
u_tip = []
v_tip = []
time.append(0.0)
u_tip.append(0.0)
v_tip.append(0.0)

while interface.is_coupling_ongoing():
    if interface.requires_writing_checkpoint():
        pass

    precice_dt = interface.get_max_time_step_size()
    dt = min([solver_dt, precice_dt])

    write_data[:, 0] = F_max * vertices[:, 1] / H  # linearly increasing load
    write_data[:, 1] = 0

    interface.write_data(mesh_name, write_data_name, vertex_ids, write_data)
    interface.advance(dt)

    if interface.requires_reading_checkpoint():
        pass

interface.finalize()
