"""
The basic example is taken from "Langtangen, Hans Petter, and Anders Logg. Solving PDEs in Python: The FEniCS
Tutorial I. Springer International Publishing, 2016."

The example code has been extended with preCICE API calls and mixed boundary conditions to allow for a Dirichlet-Neumann
coupling of two separate heat equations.

The original source code can be found on https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft03_heat.py.

Heat equation with Dirichlet conditions. (Dirichlet problem)
  u'= Laplace(u) + f  in the unit square [0,1] x [0,1]
  u = u_C             on the coupling boundary at x = 1
  u = u_D             on the remaining boundary
  u = u_0             at t = 0
  u = 1 + x^2 + alpha*y^2 + \beta*t
  f = beta - 2 - 2*alpha

Heat equation with mixed boundary conditions. (Neumann problem)
  u'= Laplace(u) + f  in the shifted unit square [1,2] x [0,1]
  du/dn = f_N         on the coupling boundary at x = 1
  u = u_D             on the remaining boundary
  u = u_0             at t = 0
  u = 1 + x^2 + alpha*y^2 + \beta*t
  f = beta - 2 - 2*alpha
"""

from __future__ import print_function
from fenics import *
import numpy as np
from enum import Enum
from matplotlib import pyplot as plt


class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem


class ComplementaryBoundary(SubDomain):
    def __init__(self, subdomain):
        self.complement = subdomain
        SubDomain.__init__(self)
    def inside(self, x, on_boundary):
        tol = 1E-14
        if on_boundary and not self.complement.inside(x, on_boundary):
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


problem = ProblemType.DIRICHLET
T = 2.0  # final time
num_steps = 10  # number of time steps
dt = T / num_steps  # time step size
alpha = 3  # parameter alpha
beta = 1.2  # parameter beta
x_coupling = 1  # x coordinate of coupling interface

# Create mesh and define function space
nx = ny = 8

if problem is ProblemType.DIRICHLET:
    p0 = Point(0, 0)
    p1 = Point(1, 1)
elif problem is ProblemType.NEUMANN:
    p0 = Point(1, 0)
    p1 = Point(2, 1)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t',
                 degree=2, alpha=alpha, beta=beta, t=0)
# Define flux on coupling interface (grad(u_D) in normal direction)
f_N = Expression('2 * x[0]', degree=1)

"""
coupling_interface = CouplingBoundary()
remaining_boundary = ComplementaryBoundary(coupling_interface))

PRECICE_INIT_MESH
PRECICE_INIT_BC
"""

def boundary(x, on_boundary):
    return on_boundary

bcs = []
bcs.append(DirichletBC(V, u_D, boundary))

# Define initial value
u_n = interpolate(u_D, V)
# u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2 * alpha)

F = u * v * dx + dt * dot(grad(u), grad(v)) * dx - (u_n + dt * f) * v * dx

"""
if problem is ProblemType.DIRICHLET:
    bcs.append(DirichletBC(V, FROM_PRECICE, coupling_interface))  # apply Dirichlet boundary condition on coupling interface
if problem is ProblemType.NEUMANN:
    F += FROM_PRECICE    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
"""

a, L = lhs(F), rhs(F)

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):
    # Update current time
    t += dt
    u_D.t = t

    PRECICE_COUPLING_NOT_CONVERGED = True

    while PRECICE_COUPLING_NOT_CONVERGED:
        """
        if problem is ProblemType.DIRICHLET:
            READ_UPDATE_BC
        if problem is ProblemType.NEUMANN:
            READ_UPDATE_BC
        """
        # Compute solution
        solve(a == L, u, bcs)
        """
        if problem is ProblemType.DIRICHLET:
            WRITE_UPDATE_BC
        if problem is ProblemType.NEUMANN:
            WRITE_UPDATE_BC
        """
        PRECICE_COUPLING_NOT_CONVERGED = False  # todo currently only one iteration

    # Plot solution
    plot(u)
    plt.pause(.1)

    # Compute error at vertices
    u_e = interpolate(u_D, V)
    error = np.abs(u_e.vector().get_local() - u.vector().get_local()).max()
    print('t = %.2f: error = %.3g' % (t, error))

    # Update previous solution
    u_n.assign(u)

# Hold plot
plt.show()
