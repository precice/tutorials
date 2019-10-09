"""
Problem setup for HT fenics-fenics tutorial
"""

from fenics import SubDomain, Point, RectangleMesh, near, Function, VectorFunctionSpace, Expression
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


def get_problem_setup(args):
    if args.dirichlet:
        problem = ProblemType.DIRICHLET
    if args.neumann:
        problem = ProblemType.NEUMANN
    if args.dirichlet and args.neumann:
        raise Exception("you can only choose either a dirichlet problem (option -d) or a neumann problem (option -n)")
    if not (args.dirichlet or args.neumann):
        raise Exception("you have to choose either a dirichlet problem (option -d) or a neumann problem (option -n)")

    if args.interface == ['simple']:
        if args.domain == ['left']:
            return DomainPart.LEFT, problem
        if args.domain == ['right']:
            return DomainPart.RIGHT, problem
        if args.dirichlet and args.neumann:
            raise Exception("you can only choose to either compute the left part of the domain (option -dl) or the right part (option -dr)")
        if not (args.domain == ['left'] or args.domain == ['right']):
            print("Default domain partitioning for simple interface is used: Left part of domain is a Dirichlet-type problem; right part is a Neumann-type problem")
            if problem is ProblemType.DIRICHLET:
                return DomainPart.LEFT, problem
            elif problem is ProblemType.NEUMANN:
                return DomainPart.RIGHT, problem

    if args.interface == ['complex']:
        if args.domain == ['circular']:
            return DomainPart.CIRCULAR, problem
        if args.domain == ['rectangle']:
            return DomainPart.RECTANGLE, problem
        if args.dirichlet and args.neumann:
            raise Exception("you can only choose to either compute the circular part of the domain (option -dc) or the residual part (option -dnc)")
        if not (args.domain == ['circular'] or args.domain == ['rectangle']):
            print("Default domain partitioning for complex interface is used: Circular part of domain is a Dirichlet-type problem; Rest of the domain is a Neumann-type problem")
            if problem is ProblemType.DIRICHLET:
                return DomainPart.CIRCULAR, problem
            elif problem is ProblemType.NEUMANN:
                return DomainPart.RECTANGLE, problem


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
    elif domain_part is DomainPart.RECTANGLE:
        n_vertices = n_vertices

    if domain_part is DomainPart.LEFT or domain_part is DomainPart.RIGHT:
        if domain_part is DomainPart.LEFT:
            p0 = Point(x_left, y_bottom)
            p1 = Point(x_coupling, y_top)
        elif domain_part is DomainPart.RIGHT:
            p0 = Point(x_coupling, y_bottom)
            p1 = Point(x_right, y_top)
        mesh = RectangleMesh(p0, p1, nx, ny)

    if domain_part is DomainPart.CIRCULAR or domain_part is DomainPart.RECTANGLE:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_right, y_top)
        whole_domain = mshr.Rectangle(p0, p1)
        if domain_part is DomainPart.CIRCULAR:
            circular_domain = mshr.Circle(midpoint, radius, n_vertices)
            mesh = mshr.generate_mesh(circular_domain, high_resolution, "cgal")
        elif domain_part is DomainPart.RECTANGLE:
            circular_domain = mshr.Circle(midpoint, radius, n_vertices)
            mesh = mshr.generate_mesh(whole_domain - circular_domain, low_resolution, "cgal")

    coupling_boundary = CouplingBoundary()
    coupling_boundary.get_user_input_args(args)

    remaining_boundary = OuterBoundary()
    remaining_boundary.get_user_input_args(args)

    return mesh, coupling_boundary, remaining_boundary


def get_facet_normal(mesh):
    """
    Manually calculate FacetNormal function
    :param mesh: fenics mesh
    :return: vector of facet normals of mesh
    """
    """
    if not mesh.type().dim() == 1:
        raise ValueError('Only works for 2-D mesh')
    """
    vertices = mesh.coordinates()
    cells = mesh.cells()

    vec1 = vertices[cells[:, 1]] - vertices[cells[:, 0]]
    normals = vec1[:, [1, 0]] * np.array([1, -1])
    normals /= np.sqrt((normals**2).sum(axis=1))[:, np.newaxis]

    # Ensure outward pointing normal
    mesh.init_cell_orientations(Expression(('x[0]', 'x[1]'), degree=1))
    normals[mesh.cell_orientations() == 1] *= -1

    V = VectorFunctionSpace(mesh, 'DG', 0)
    norm = Function(V)
    nv = norm.vector()

    for n in (0, 1):
        dofmap = V.sub(n).dofmap()
        for i in range(dofmap.global_dimension()):
            dof_indices = dofmap.cell_dofs(i)
            assert len(dof_indices) == 1
            nv[dof_indices[0]] = normals[i, n]

    return norm
