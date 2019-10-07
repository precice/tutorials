"""
Problem setup for HT fenics-fenics tutorial
"""

from fenics import SubDomain, Point, RectangleMesh, near
from my_enums import DomainPart, ProblemType
import mshr

y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface
radius = 0.2
midpoint = Point(0.5, 0.5)


class OuterBoundary(SubDomain):
    def get_user_input_args(self, args):
        self._interface = args.interface

    def inside(self, x, on_boundary):
        tol = 1E-14
        if self._interface == ['simple']:
            if on_boundary and not near(x[0], x_coupling, tol) or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
                return True
            else:
                return False
        if self._interface == ['complex']:
            point = Point(x[0], x[1])
            if on_boundary and point.distance(midpoint)**2 - radius**2 > tol or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
                return True
            else:
                return False


class CouplingBoundary(SubDomain):
    def get_user_input_args(self, args):
        self._interface = args.interface

    def inside(self, x, on_boundary):
        tol = 1E-14
        if self._interface == ['simple']:
            if on_boundary and near(x[0], x_coupling, tol):
                return True
            else:
                return False
        if self._interface == ['complex']:
            point = Point(x[0], x[1])
            if on_boundary and point.distance(midpoint)**2 - radius**2 < tol:
                return True
            else:
                return False


def get_problem_setup(args, problem):
    if args.interface == ['simple']:
        if args.domain == ['left']:
            print("Reaching command return DomainPart.LEFT")
            return DomainPart.LEFT
        if args.domain == ['right']:
            return DomainPart.RIGHT
        if args.dirichlet and args.neumann:
            raise Exception("you can only choose to either compute the left part of the domain (option -dl) or the right part (option -dr)")
        if not (args.domain == ['left'] or args.domain == ['right']):
            print("Default domain partitioning for simple interface is used: Left part of domain is a Dirichlet-type problem; right part is a Neumann-type problem")
            if problem is ProblemType.DIRICHLET:
                return DomainPart.LEFT
            elif problem is ProblemType.NEUMANN:
                return DomainPart.RIGHT

    if args.interface == ['complex']:
        if args.domain == ['circular']:
            return DomainPart.CIRCULAR
        if args.domain == ['rest']:
            return DomainPart.REST
        if args.dirichlet and args.neumann:
            raise Exception("you can only choose to either compute the circular part of the domain (option -dc) or the residual part (option -dnc)")
        if not (args.domain == ['circular'] or args.domain == ['rest']):
            print("Default domain partitioning for complex interface is used: Circular part of domain is a Dirichlet-type problem; Rest of the domain is a Neumann-type problem")
            if problem is ProblemType.DIRICHLET:
                return DomainPart.CIRCULAR
            elif problem is ProblemType.NEUMANN:
                return DomainPart.REST


def get_geometry(args, domain_part):
    nx = 5
    ny = 10
    low_resolution = 5
    high_resolution = 5
    n_vertices = 3

    if domain_part is DomainPart.LEFT:
        nx = nx*3
    elif domain_part is DomainPart.RIGHT:
        ny = 20

    if domain_part is DomainPart.CIRCULAR:
        n_vertices = n_vertices
    elif domain_part is DomainPart.REST:
        n_vertices = n_vertices

    if args.interface == ['simple']:
        if domain_part is DomainPart.LEFT:
            p0 = Point(x_left, y_bottom)
            p1 = Point(x_coupling, y_top)
        elif domain_part is DomainPart.RIGHT:
            p0 = Point(x_coupling, y_bottom)
            p1 = Point(x_right, y_top)
        return RectangleMesh(p0, p1, nx, ny)

    if args.interface == ['complex']:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_right, y_top)
        whole_domain = mshr.Rectangle(p0, p1)
        if domain_part is DomainPart.CIRCULAR:
            circular_domain = mshr.Circle(midpoint, radius, n_vertices)
            return mshr.generate_mesh(circular_domain, high_resolution, "cgal")
        elif domain_part is DomainPart.REST:
            circular_domain = mshr.Circle(midpoint, radius, n_vertices)
            return mshr.generate_mesh(whole_domain - circular_domain, low_resolution, "cgal")

