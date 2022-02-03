"""
Problem setup for HT fenics-fenics tutorial
"""

from fenics import SubDomain, Point, RectangleMesh, near, Function, VectorFunctionSpace, Expression
from my_enums import DomainPart, ProblemType
import mshr
import numpy as np

y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface
radius = 0.2
midpoint = Point(0.5, 0.5)


class ExcludeStraightBoundary(SubDomain):
    def get_user_input_args(self, args):
        self._interface = args.interface

    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not near(x[0], x_coupling, tol) or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
            return True
        else:
            return False


class ExcludeCircleBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        point = Point(x[0], x[1])
        if on_boundary and not point.distance(midpoint)**2 - radius**2 < tol:
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


class CircleBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        point = Point(x[0], x[1])
        if on_boundary and point.distance(midpoint)**2 - radius**2 < tol:
            return True
        else:
            return False


def get_problem_setup(args):
    if args.dirichlet and not args.neumann:
        problem = ProblemType.DIRICHLET
    elif args.neumann and not args.dirichlet:
        problem = ProblemType.NEUMANN
    elif args.dirichlet and args.neumann:
        raise Exception("you can only choose either a dirichlet problem (option -d) or a neumann problem (option -n)")
    elif not (args.dirichlet or args.neumann):
        raise Exception("you have to choose either a dirichlet problem (option -d) or a neumann problem (option -n)")
    else:
        raise Exception("invalid control flow.")

    if args.interface == 'simple':
        if args.domain == 'circular' or args.domain == 'rectangle':
            raise Exception("Only --domain left or --domain right is supported for --interface {}. "
                            "Invalid --domain {} provided.".format(args.interface, args.domain))
        elif args.domain == 'left':
            return DomainPart.LEFT, problem
        elif args.domain == 'right':
            return DomainPart.RIGHT, problem
        elif args.dirichlet and args.neumann:
            raise Exception("You can only choose to either compute the left part of the domain (option -dl) "
                            "or the right part (option -dr)")
        elif not (args.domain == 'left' or args.domain == 'right'):
            print("Default domain partitioning for simple interface is used: Left part of domain is a "
                  "Dirichlet-type problem; right part is a Neumann-type problem")
            if problem is ProblemType.DIRICHLET:
                return DomainPart.LEFT, problem
            elif problem is ProblemType.NEUMANN:
                return DomainPart.RIGHT, problem
    elif args.interface == 'complex':
        if args.domain == 'left' or args.domain == 'right':
            raise Exception("Only --domain circular or --domain rectangle is supported for --interface {}. "
                            "Invalid --domain {} provided.".format(args.interface, args.domain))
        elif args.domain == 'circular':
            return DomainPart.CIRCULAR, problem
        elif args.domain == 'rectangle':
            return DomainPart.RECTANGLE, problem
        elif args.dirichlet and args.neumann:
            raise Exception("You can only choose to either compute the circular part of the domain (option -dc) "
                            "or the residual part (option -dnc)")
        elif not (args.domain == 'circular' or args.domain == 'rectangle'):
            print("Default domain partitioning for complex interface is used: Circular part of domain "
                  "is a Neumann-type problem; Rest of the domain is a Dirichlet-type problem")
            if problem is ProblemType.NEUMANN:
                return DomainPart.CIRCULAR, problem
            elif problem is ProblemType.DIRICHLET:
                return DomainPart.RECTANGLE, problem
    else:
        raise Exception("invalid interface provided: args.interface = {}".format(args.interface))


def get_geometry(domain_part):
    nx = 5
    ny = 10
    low_resolution = 5
    high_resolution = 5
    n_vertices = 20

    if domain_part is DomainPart.LEFT:
        nx = nx * 3
    elif domain_part is DomainPart.RIGHT:
        ny = ny * 2
    elif domain_part is DomainPart.CIRCULAR:
        n_vertices = n_vertices
    elif domain_part is DomainPart.RECTANGLE:
        n_vertices = n_vertices
    else:
        raise Exception("invalid domain_part: {}".format(domain_part))

    if domain_part is DomainPart.LEFT or domain_part is DomainPart.RIGHT:
        if domain_part is DomainPart.LEFT:
            p0 = Point(x_left, y_bottom)
            p1 = Point(x_coupling, y_top)
        elif domain_part is DomainPart.RIGHT:
            p0 = Point(x_coupling, y_bottom)
            p1 = Point(x_right, y_top)
        else:
            raise Exception("invalid control flow!")
        mesh = RectangleMesh(p0, p1, nx, ny)
        coupling_boundary = StraightBoundary()
        remaining_boundary = ExcludeStraightBoundary()

    elif domain_part is DomainPart.CIRCULAR or domain_part is DomainPart.RECTANGLE:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_right, y_top)
        whole_domain = mshr.Rectangle(p0, p1)
        if domain_part is DomainPart.CIRCULAR:
            circular_domain = mshr.Circle(midpoint, radius, n_vertices)
            mesh = mshr.generate_mesh(circular_domain, high_resolution, "cgal")
        elif domain_part is DomainPart.RECTANGLE:
            circular_domain = mshr.Circle(midpoint, radius, n_vertices)
            mesh = mshr.generate_mesh(whole_domain - circular_domain, low_resolution, "cgal")
        else:
            raise Exception("invalid control flow!")
        coupling_boundary = CircleBoundary()
        remaining_boundary = ExcludeCircleBoundary()

    else:
        raise Exception("invalid control flow!")

    return mesh, coupling_boundary, remaining_boundary
