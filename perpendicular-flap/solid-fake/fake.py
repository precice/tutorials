from __future__ import division

import numpy as np
import precice


def displace_flap(x, y, t, flap_tip_y):
    x_displ = np.zeros_like(x)
    y_displ = np.zeros_like(y)
    # first, get displacement independent of x, only dependent on y and t
    # maximal displacement at flap tip should be 0.5
    # initially, flap's displacement is 0
    x_displ = np.minimum(((np.sin(3 * np.pi * t + np.arcsin(-0.95)) + 0.95) / 8) * y / flap_tip_y, 0.5 * y / flap_tip_y)

    displ = np.zeros((len(x), 2))
    displ[:, 0] = x_displ
    displ[:, 1] = y_displ

    return displ


configuration_file_name = "../precice-config.xml"
participant_name = "Solid"
mesh_name = "Solid-Mesh"
write_data_name = 'Displacement'

solver_process_index = 0
solver_process_size = 1

# define mesh
H = 1
W = 0.1

interface = precice.Participant(participant_name, configuration_file_name, solver_process_index, solver_process_size)
dimensions = interface.get_mesh_dimensions(mesh_name)
assert (dimensions == 2)

x_left = 0.0 - 0.5 * W  # left boundary of the flap
x_right = 0.5 * W  # right boundary of the flap
y_bottom = 0.0  # bottom of the flap
y_top = y_bottom + H  # top of the flap

n = 24  # Number of vertices per side
t = 0

vertices = np.zeros((2 * n, dimensions))
vertices[:n, 1] = np.linspace(y_bottom, y_top, n)
vertices[n:, 1] = np.linspace(y_bottom, y_top, n)
vertices[:n, 0] = x_left
vertices[n:, 0] = x_right

vertex_ids = interface.set_mesh_vertices(mesh_name, vertices)

interface.initialize()
# change if necessary
solver_dt = np.inf
# for checkpointing
t_cp = 0

while interface.is_coupling_ongoing():
    if interface.requires_writing_checkpoint():
        t_cp = t

    precice_dt = interface.get_max_time_step_size()
    dt = min([solver_dt, precice_dt])
    # wiggle the flap
    write_data = displace_flap(vertices[:, 0], vertices[:, 1], t, H)

    interface.write_data(mesh_name, write_data_name, vertex_ids, write_data)
    interface.advance(dt)

    if interface.requires_reading_checkpoint():
        t = t_cp
    else:
        # update t
        t += dt

interface.finalize()
