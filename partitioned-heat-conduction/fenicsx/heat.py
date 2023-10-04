"""
The basic example is taken from "Langtangen, Hans Petter, and Anders Logg. Solving PDEs in Python: The FEniCS
Tutorial I. Springer International Publishing, 2016."

The example code has been extended with preCICE API calls and mixed boundary conditions to allow for a Dirichlet-Neumann
coupling of two separate heat equations. It also has been adapted to be compatible with FEniCSx.

The original source code can be found on https://jsdokken.com/dolfinx-tutorial/chapter2/heat_code.html.

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

from __future__ import print_function, division
from mpi4py import MPI
from dolfinx.fem import Function, FunctionSpace, Expression, Constant, dirichletbc, locate_dofs_geometrical
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile
from ufl import TrialFunction, TestFunction, dx, ds, dot, grad, inner, lhs, rhs, FiniteElement, VectorElement
from fenicsxprecice import Adapter
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
import argparse
import numpy as np
from problem_setup import get_geometry


def determine_gradient(V_g, u):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(grad(u), v) * dx
    problem = LinearProblem(a, L)
    return problem.solve()


parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
command_group = parser.add_mutually_exclusive_group(required=True)
command_group.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest="dirichlet",
                           action="store_true")
command_group.add_argument("-n", "--neumann", help="create a neumann problem", dest="neumann", action="store_true")
parser.add_argument("-e", "--error-tol", help="set error tolerance", type=float, default=10**-6,)

args = parser.parse_args()

fenics_dt = .1  # time step size
# Error is bounded by coupling accuracy. In theory we would obtain the analytical solution.
error_tol = args.error_tol

alpha = 3  # parameter alpha
beta = 1.2  # parameter beta

if args.dirichlet and not args.neumann:
    problem = ProblemType.DIRICHLET
    domain_part = DomainPart.LEFT
elif args.neumann and not args.dirichlet:
    problem = ProblemType.NEUMANN
    domain_part = DomainPart.RIGHT

mesh, coupling_boundary, remaining_boundary = get_geometry(MPI.COMM_WORLD, domain_part)

# Define function space using mesh
scalar_element = FiniteElement("P", mesh.ufl_cell(), 2)
vector_element = VectorElement("P", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, scalar_element)
V_g = FunctionSpace(mesh, vector_element)
W = V_g.sub(0).collapse()[0]

# Define boundary conditions


class Expression_u_D:
    def __init__(self):
        self.t = 0.0

    def eval(self, x):
        return np.full(x.shape[1], 1 + x[0] * x[0] + alpha * x[1] * x[1] + beta * self.t)


u_D = Expression_u_D()
u_D_function = Function(V)
u_D_function.interpolate(u_D.eval)

if problem is ProblemType.DIRICHLET:
    # Define flux in x direction
    class Expression_f_N:
        def __init__(self):
            pass

        def eval(self, x):
            return np.full(x.shape[1], 2 * x[0])

    f_N = Expression_f_N()
    f_N_function = Function(V)
    f_N_function.interpolate(f_N.eval)

# Define initial value
u_n = Function(V, name="Temperature")
u_n.interpolate(u_D.eval)

precice, precice_dt, initial_data = None, 0.0, None

# Initialize the adapter according to the specific participant
if problem is ProblemType.DIRICHLET:
    precice = Adapter(MPI.COMM_WORLD, adapter_config_filename="precice-adapter-config-D.json")
    precice_dt = precice.initialize(coupling_boundary, read_function_space=V, write_object=f_N_function)
elif problem is ProblemType.NEUMANN:
    precice = Adapter(MPI.COMM_WORLD, adapter_config_filename="precice-adapter-config-N.json")
    precice_dt = precice.initialize(coupling_boundary, read_function_space=W, write_object=u_D_function)

dt = Constant(mesh, 0.0)
dt.value = np.min([fenics_dt, precice_dt])

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)


class Expression_f:
    def __init__(self):
        self.t = 0.0

    def eval(self, x):
        return np.full(x.shape[1], beta - 2 - 2 * alpha)


f = Expression_f()
f_function = Function(V)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f_function) * v * dx

dofs_remaining = locate_dofs_geometrical(V, remaining_boundary)
bcs = [dirichletbc(u_D_function, dofs_remaining)]

coupling_expression = precice.create_coupling_expression()
read_data = precice.read_data()
precice.update_coupling_expression(coupling_expression, read_data)
if problem is ProblemType.DIRICHLET:
    # modify Dirichlet boundary condition on coupling interface
    dofs_coupling = locate_dofs_geometrical(V, coupling_boundary)
    bcs.append(dirichletbc(coupling_expression, dofs_coupling))
if problem is ProblemType.NEUMANN:
    # modify Neumann boundary condition on coupling interface, modify weak
    # form correspondingly
    F += v * coupling_expression * ds

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V, name="Temperature")
t = 0

# reference solution at t=0
u_ref = Function(V, name="reference")
u_ref.interpolate(u_D_function)

# Generating output files
temperature_out = XDMFFile(MPI.COMM_WORLD, f"./output/{precice.get_participant_name()}.xdmf", "w")
temperature_out.write_mesh(mesh)
ref_out = XDMFFile(MPI.COMM_WORLD, f"./output/ref{precice.get_participant_name()}.xdmf", "w")
ref_out.write_mesh(mesh)

# output solution at t=0, n=0
n = 0

temperature_out.write_function(u_n, t)
ref_out.write_function(u_ref, t)

u_D.t = t + dt.value
u_D_function.interpolate(u_D.eval)
f.t = t + dt.value
f_function.interpolate(f.eval)

if problem is ProblemType.DIRICHLET:
    flux = Function(V_g, name="Heat-Flux")

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.is_action_required(precice.action_write_iteration_checkpoint()):
        precice.store_checkpoint(u_n, t, n)

    read_data = precice.read_data()

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(coupling_expression, read_data)

    dt.value = np.min([fenics_dt, precice_dt])

    linear_problem = LinearProblem(a, L, bcs=bcs)
    u_np1 = linear_problem.solve()

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        flux = determine_gradient(V_g, u_np1)
        flux_x = Function(W)
        flux_x.interpolate(flux.sub(0))
        precice.write_data(flux_x)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write_data(u_np1)

    precice_dt = precice.advance(dt.value)

    # roll back to checkpoint
    if precice.is_action_required(precice.action_read_iteration_checkpoint()):
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.interpolate(u_cp)
        t = t_cp
        n = n_cp
    else:  # update solution
        u_n.interpolate(u_np1)
        t += dt.value
        n += 1

    if precice.is_time_window_complete():
        u_ref.interpolate(u_D_function)
        error = compute_errors(u_n, u_ref, total_error_tol=error_tol)
        print('n = %d, t = %.2f: L2 error on domain = %.3g' % (n, t, error))
        print('output u^%d and u_ref^%d' % (n, n))

        temperature_out.write_function(u_n, t)
        ref_out.write_function(u_ref, t)

    # Update Dirichlet BC
    u_D.t = t + dt.value
    u_D_function.interpolate(u_D.eval)
    f.t = t + dt.value
    f_function.interpolate(f.eval)


# Hold plot
precice.finalize()
