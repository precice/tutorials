from __future__ import division

import numpy as np
import precice


def displace_flap(x, y, t, flap_tip_y):
    x_displ = np.zeros_like(x)
    y_displ = np.zeros_like(y)
    # get displacement independent of x, only dependent on y and t
    max_x_displ = 0.5
    period_fac = 3 * np.pi
    damping_fac = 8  # damps the amplitude of the sine
    # defines how much the sine is shifted in y-direction
    shift = 0.95
    # wiggles the flap periodically.
    # the arcsin(-shift) in the sine evaluation is necessary to start at a flap displacement of 0 at t=0
    # (i.e. sin(arcsin(-shift))+shift = 0)
    x_displ = np.minimum(((np.sin(period_fac * t + np.arcsin(-shift)) + shift) /
                         damping_fac), max_x_displ) * y / flap_tip_y

    dimensions = 2
    displ = np.zeros((len(x), dimensions))
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
# define vertices of flap's left side
vertices[:n, 1] = np.linspace(y_bottom, y_top, n)
vertices[:n, 0] = x_left
# define vertices of flap's right side
vertices[n:, 1] = np.linspace(y_bottom, y_top, n)
vertices[n:, 0] = x_right

vertex_ids = interface.set_mesh_vertices(mesh_name, vertices)

if interface.requires_initial_data():
    # initially, there should be no displacement
    interface.write_data(np.zeros_like(vertices))

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
