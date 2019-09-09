import sympy as sp
from my_enums import DomainPart, ProblemType
from fenics import Expression, RectangleMesh, Point, SubDomain, near

y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface


class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not near(x[0], x_coupling, tol) or near(x[1], y_top, tol) or near(x[1], y_bottom, tol):
            return True
        else:
            return False


class CouplingBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and near(x[0], x_coupling, tol):
            return True
        else:
            return False


def get_problem_setup(args, domain_part, problem):

    alpha = args.alpha  # parameter alpha
    beta = args.beta  # parameter beta
    gamma = args.gamma  # parameter gamma, dependence of heat flux on time

    # Define boundary condition
    x, y, t, dt = sp.symbols('x[0], x[1], t, dt')

    # Define analytical solution
    if args.time_dependence == "l":
        p = 1
        u_analytical = 1 + gamma * (1 + t)**p * x**2 + (1-gamma) * x**2 + alpha * y**2 + beta * t
    if args.time_dependence == "q":
        p = 2
        u_analytical = 1 + gamma * (1 + t)**p * x**2 + (1-gamma) * x**2 + alpha * y**2 + beta * t
    elif args.time_dependence == "s":
        u_analytical = 1 + gamma * sp.sin(t) * x**2 + (1-gamma) * x**2 + alpha * y**2 + beta * t
    u_D = Expression(sp.printing.ccode(u_analytical), degree=2, t=0)

    # Compute corresponding right-hand side
    f_rhs = - sp.diff(sp.diff(u_analytical, x), x) - sp.diff(sp.diff(u_analytical, y), y) + sp.diff(u_analytical, t)
    f_np1 = Expression(sp.printing.ccode(f_rhs), degree=2, t=0)
    f_n = Expression(sp.printing.ccode(f_rhs), degree=2, t=0)

    # Compute flux in x direction on coupling interface (grad(u_D) in normal direction)
    flux_analytical = sp.diff(u_analytical, x)
    if (domain_part is DomainPart.LEFT and problem is ProblemType.DIRICHLET) or \
            (domain_part is DomainPart.RIGHT and problem is ProblemType.NEUMANN):
        f_N = Expression(sp.printing.ccode(flux_analytical), degree=1, t=0)
    elif (domain_part is DomainPart.RIGHT and problem is ProblemType.DIRICHLET) or \
            (domain_part is DomainPart.LEFT and problem is ProblemType.NEUMANN):
        f_N = Expression(sp.printing.ccode(-flux_analytical), degree=1, t=0)

    return f_np1, f_n, u_D, f_N


def get_geometry(domain_part):

    nx = 10
    ny = 10

    if domain_part is DomainPart.LEFT:
        nx = nx * 3
    elif domain_part is DomainPart.RIGHT:
        ny = 20

    if domain_part is DomainPart.LEFT:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_coupling, y_top)
    elif domain_part is DomainPart.RIGHT:
        p0 = Point(x_coupling, y_bottom)
        p1 = Point(x_right, y_top)

    return RectangleMesh(p0, p1, nx, ny)
