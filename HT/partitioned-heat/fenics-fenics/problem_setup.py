import sympy as sp
from my_enums import DomainPart, ProblemType
from fenics import Expression, RectangleMesh, Point, SubDomain, near
import argparse

y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1  # x coordinate of coupling interface


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


class CompleteBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary:
            return True
        else:
            return False


def get_manufactured_solution(time_dependence, alpha, beta, gamma):
    x, y, t, dt = sp.symbols('x[0], x[1], t, dt')
    # Define analytical solution
    if time_dependence == "l":  # linear
        p = 1
        g = (t + 1)**p  # g_pol
    elif time_dependence == "q":  # quadratic
        p = 2
        g = (t + 1)**p  # g_pol
    elif time_dependence == "c":  # cubic
        p = 3
        g = (t + 1)**p  # g_pol
    elif time_dependence == "qrt":  # quartic
        p = 4
        g = (t + 1)**p  # g_pol
    elif time_dependence == "s":  # sinusoidal
        g = sp.sin(t) # g_tri
    manufactured_solution = 1 + gamma * g * x**2 + (1-gamma) * x**2 + alpha * y**2 + beta * t
    print("manufactured solution = {}".format(manufactured_solution))
    return manufactured_solution


def get_problem_setup(args, domain_part, problem):

    # Define boundary condition
    x, y, t, dt = sp.symbols('x[0], x[1], t, dt')
    u_analytical = get_manufactured_solution(args.time_dependence, args.alpha, args.beta, args.gamma)
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
    else:
        f_N = Expression("0", degree=0)

    return f_np1, f_n, u_D, f_N


def get_geometry(domain_part, nx, ny):

    if domain_part is DomainPart.LEFT:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_coupling, y_top)
    elif domain_part is DomainPart.RIGHT:
        p0 = Point(x_coupling, y_bottom)
        p1 = Point(x_right, y_top)
    elif domain_part is DomainPart.ALL:
        p0 = Point(x_left, y_bottom)
        p1 = Point(x_right, y_top)

    return RectangleMesh(p0, p1, nx, ny)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-gamma", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=1.0, type=float)
    parser.add_argument("-alpha", "--alpha", help="parameter gamma to set temporal dependence of heat flux", default=3.0, type=float)
    parser.add_argument("-beta", "--beta", help="parameter gamma to set temporal dependence of heat flux", default=1.2, type=float)
    parser.add_argument("-t", "--time-dependence", help="choose whether there is a linear (l), quadratic (q) or sinusoidal (s) dependence on time", type=str, default="l")
    parser.add_argument("-mth", "--method", help="time stepping method being used", default='ie')
    args = parser.parse_args()
    u = get_manufactured_solution(args.time_dependence, args.alpha, args.beta, args.gamma)
    u = u.subs('x[0]', 'x').subs('x[1]','y')
    print(sp.latex(u))

