"""
Problem setup for partitioned-heat-conduction/fenics tutorial
"""

import sympy as sp
from fenics import SubDomain, Point, RectangleMesh, near, Function, Expression
from my_enums import DomainPart, TimeDependence


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
    nx = ny = 15

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


def get_manufactured_solution(time_dependence, alpha, beta, p=-1):
    x, y, t = sp.symbols('x[0], x[1], t')
    # Define analytical solution
    if time_dependence == TimeDependence.POLYNOMIAL:
        assert (p > -1)
        g = (t + 1)**p
    elif time_dependence == TimeDependence.TRIGONOMETRIC:
        g = 1 + sp.sin(t)
    else:
        raise Exception(f"Invalid TimeDependence {time_dependence} provided.")

    manufactured_solution = 1 + g * x**2 + alpha * y**2 + beta * t
    print("manufactured solution = {}".format(manufactured_solution))
    return manufactured_solution, {'x': x, 'y': y, 't': t}
