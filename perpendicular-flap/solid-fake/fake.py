from __future__ import division

import numpy as np
import precice
def displacement_param(t):
    # this parameter exists to define a time-dependent displacement of the flap's tip
    if t > 1:
        return 0.5
    else:
        return 0.5*t


# computes the flap displacements
def compute_flap_displacements(t, height_flap, width_flap, vertices_per_side, dimensions=2):
    a = displacement_param(t)
    # ASSUMPTION: the middle of the tip must still have a length of 1
    # For that, one can use the formula height_flap=int_0^y(sqrt(1+f'(z)^2) dz
    #
    # In the following, the function 1/a cosh(ay) - 1/a is used because the analytical solution can be determined for
    # the height_flap
    #
    # We get for the y-coordinate of the tip's middle: y=arcsinh(a*height_flap)/a
    # the y-displacement is therefore: y_displ =height_flap - y
    # the x-coordinates can hence be computed by x(y)= 1/a cosh(ay) - 1/a in [0, y]

    # for a = 0, we just use a straight line
    if a == 0:
        # if a=0 the flap is straight -> displacements are all 0
        return np.zeros((2 * vertices_per_side, dimensions))

    fun_x_coord = lambda y: (np.cosh(a * y) - 1) / a

    ##########################
    # IN THE FOLLOWING: a!=0 #
    ##########################

    #uniformly distributed vertices
    h_vert = height_flap / (vertices_per_side - 1)
    # with the function above, we know, that the middle flap's length can be computed with len=sinh(a*x)/a
    # now, we will compute the points of the flap_shape, were the len is not the height_flap, but n*h_vert
    # the reason: we guarantee that we compute the points of the borders uniformly distributed on the basis that
    # the flap is straight => n*h_vert = sinh(a*x)/a <=> x = arcsinh(n*h_vert*a)/a
    displacements_right = np.zeros((vertices_per_side, dimensions))
    displacements_left = np.zeros((vertices_per_side, dimensions))

    # get uniformly distanced (x,y) points of the flap
    y_values_middle = np.array([np.arcsinh(i * h_vert * a) / a for i in range(vertices_per_side)])
    x_values_middle = np.array([fun_x_coord(y_values_middle[i]) for i in range(len(y_values_middle))])
    # dummy values for [0] entries to avoid warnings
    x_values_middle[0] = y_values_middle[0] = 1

    # now we compute the (x,y) values of the border according to our (x,y) values of the flap's middle

    # step 1: compute the tangent with slope normal to the flap's middle shape
    # derivative of fun_x_coord is sinh(ay)
    slopes = -1 / np.sinh(a * y_values_middle)

    # step 2: compute the x and y coordinates that are on the tangent of step 1 and distanced width_flap/2 from
    # (x_values_middle,y_values_middle)
    # use slope = dx/dy and (width/2)^2 = dy^2 + dx^2
    # -> dy^2 = (width/2)^2/(slope^2+1)
    # we get two results from that because the result is the absolute value!
    dy = np.divide(width_flap / 2, np.sqrt(slopes * slopes + 1))
    # -> dx = slope*dy (second dx only negative value of result)
    dx = slopes * dy

    # step 3: now, we can compute with dx, dy, x_values_middle and y_values_middle the coordinates
    # of the border's vertices and eventually the displacements!
    displacements_right[:, 0] = x_values_middle - dx
    displacements_right[:, 1] = y_values_middle - dy
    displacements_left[:, 0] = x_values_middle + dx  # second dx has inverted sign
    displacements_left[:, 1] = y_values_middle + dy

    # displacements: remember that y-coordinates of the vertices are normally uniformly distributed,
    # i.e. y_0=0, y_1=h_vert, y_2=2*h_vert,...
    # and the x-coordinates are for the left border -width/2 and for the right border width/2
    default_x_left = -width_flap / 2
    default_x_right = width_flap / 2
    for i in range(vertices_per_side):
        #displacements right
        displacements_right[i][0] = displacements_right[i][0] - default_x_right
        displacements_right[i][1] = displacements_right[i][1] - h_vert * i
        #displacements left
        displacements_left[i][0] = displacements_left[i][0] - default_x_left
        displacements_left[i][1] = displacements_left[i][1] - h_vert * i
    # bottom is fixed -> x & y displacement 0
    displacements_left[0][:] = 0
    displacements_right[0][:] = 0

    return np.concatenate((displacements_left, displacements_right), axis=0)

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
vertices[:n, 0] = x_left  # vertices at left side of flap
vertices[n:, 0] = x_right  # vertices at right side of flap
vertices[:n, 1] = np.linspace(y_bottom, y_top, n)  # equally disrtibuted vertices left
vertices[n:, 1] = np.linspace(y_bottom, y_top, n)  # equally disrtibuted vertices right
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
    write_data = compute_flap_displacements(t, H, W, n)

    interface.write_data(mesh_name, write_data_name, vertex_ids, write_data)
    interface.advance(dt)

    if interface.requires_reading_checkpoint():
        t = t_cp
    else:
        # update t
        t += dt

interface.finalize()

