"""
Problem setup for partitioned-heat-conduction/fenics tutorial
"""

from fenics import SubDomain, Point, RectangleMesh, near, Function, Expression
from my_enums import DomainPart


y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.0  # x coordinate of coupling interface


class ExcludeStraightBoundary(SubDomain):
    def get_user_input_args(self, args):
        self._interface = args.interface

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not near(x[0], x_coupling, tol) or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
            return True
        else:
            return False


class StraightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[0], x_coupling, tol):
            return True
        else:
            return False


def get_geometry(domain_part):
    nx = ny = 11

    if domain_part is DomainPart.LEFT:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_coupling, y_top)
    elif domain_part is DomainPart.RIGHT:
        p0 = Point(x_coupling, y_bottom)
        p1 = Point(x_right, y_top)
    else:
        raise Exception("invalid domain_part: {}".format(domain_part))

    mesh = RectangleMesh(p0, p1, nx, ny, diagonal="left")
    coupling_boundary = StraightBoundary()
    remaining_boundary = ExcludeStraightBoundary()

    return mesh, coupling_boundary, remaining_boundary
