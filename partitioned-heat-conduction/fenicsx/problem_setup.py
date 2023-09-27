"""
Problem setup for partitioned-heat-conduction/fenicsx tutorial
"""
from dolfinx.mesh import DiagonalType, create_rectangle
from my_enums import DomainPart
import numpy as np


y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.0  # x coordinate of coupling interface


def exclude_straight_boundary(x):
    tol = 1E-14
    return np.logical_or(
        np.logical_or(~np.isclose(x[0], x_coupling, tol), np.isclose(x[1], y_top, tol)),
        np.isclose(x[1], y_bottom, tol)
    )


def straight_boundary(x):
    tol = 1E-14
    return np.isclose(x[0], x_coupling, tol)


def get_geometry(mpi_comm, domain_part):
    nx = ny = 9
    low_resolution = 5
    high_resolution = 5
    n_vertices = 20

    if domain_part is DomainPart.LEFT:
        p0 = (x_left, y_bottom)
        p1 = (x_coupling, y_top)
    elif domain_part is DomainPart.RIGHT:
        p0 = (x_coupling, y_bottom)
        p1 = (x_right, y_top)
    else:
        raise Exception("invalid domain_part: {}".format(domain_part))

    mesh = create_rectangle(mpi_comm, [p0, p1], [nx, ny], diagonal=DiagonalType.left)
    coupling_boundary = straight_boundary
    remaining_boundary = exclude_straight_boundary

    return mesh, coupling_boundary, remaining_boundary
