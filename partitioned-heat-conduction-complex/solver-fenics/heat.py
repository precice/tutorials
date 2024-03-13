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

from __future__ import print_function, division
from fenics import Function, FunctionSpace, Expression, Constant, DirichletBC, TrialFunction, TestFunction, \
    File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace, FacetNormal, MeshFunction, MPI
from fenicsprecice import Adapter
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
import argparse
import numpy as np
from problem_setup import get_geometry, get_problem_setup
import sympy as sp


def determine_gradient(V_g, u, flux):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    :param flux: returns calculated flux into this value
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(grad(u), v) * dx
    solve(a == L, flux)


parser = argparse.ArgumentParser(
    description='Solving heat equation for simple or complex interface case')
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')
parser.add_argument("-g", "--gamma", help="parameter gamma to set temporal dependence of heat flux", default=0.0,
                    type=float)
parser.add_argument("-i", "--interface", metavar="interface_type string", type=str, choices=['simple', 'complex'],
                    help="Type of coupling interface case to be solved. Options: simple, complex", default="simple")
parser.add_argument("-dom", "--domain", metavar='domain_type string', type=str,
                    choices=['left', 'right', 'circular', 'rectangle'],
                    help="Specifying part of the domain being solved. "
                    "For simple interface the options are left, right, "
                    "for complex interface the options are circular, rest")
args = parser.parse_args()

fenics_dt = .1
# Error is bounded by coupling accuracy. In theory we can obtain the analytical solution.
error_tol = 10 ** -6
alpha = 3  # parameter alpha
beta = 1.2  # parameter beta
gamma = args.gamma  # parameter gamma, dependence of heat flux on time

# Create mesh and separate mesh components for grid, boundary and coupling interface
domain_part, problem = get_problem_setup(args)
mesh, coupling_boundary, remaining_boundary = get_geometry(domain_part)

# Define function space using mesh
V = FunctionSpace(mesh, 'P', 2)
V_g = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary conditions
# create sympy expression of manufactured solution
x_sp, y_sp, t_sp = sp.symbols(['x[0]', 'x[1]', 't'])
u_D_sp = 1 + gamma * t_sp * x_sp * x_sp + (1 - gamma) * x_sp * x_sp + alpha * y_sp * y_sp + beta * t_sp
u_D = Expression(sp.ccode(u_D_sp), degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
f_N = Expression((sp.ccode(u_D_sp.diff(x_sp)), sp.ccode(u_D_sp.diff(y_sp))), degree=1, gamma=gamma, alpha=alpha, t=0)
f_N_function = interpolate(f_N, V_g)

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

precice, precice_dt, initial_data = None, 0.0, None

# Initialize the adapter according to the specific participant
precice = Adapter(adapter_config_filename="precice-adapter-config.json")

if problem is ProblemType.DIRICHLET:
    precice.initialize(coupling_boundary, read_function_space=V, write_object=f_N_function)
elif problem is ProblemType.NEUMANN:
    precice.initialize(coupling_boundary, read_function_space=V_g, write_object=u_D_function)

boundary_marker = False

dt = Constant(0)
precice_dt = precice.get_max_time_step_size()
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
# du_dt-Laplace(u) = f
f_sp = u_D_sp.diff(t_sp) - u_D_sp.diff(x_sp).diff(x_sp) - u_D_sp.diff(y_sp).diff(y_sp)
f = Expression(sp.ccode(f_sp), degree=2, alpha=alpha, beta=beta, gamma=gamma, t=0)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx

bcs = [DirichletBC(V, u_D, remaining_boundary)]

# Set boundary conditions at coupling interface once wrt to the coupling expression
coupling_expression = precice.create_coupling_expression()
if problem is ProblemType.DIRICHLET:
    # modify Dirichlet boundary condition on coupling interface
    bcs.append(DirichletBC(V, coupling_expression, coupling_boundary))
if problem is ProblemType.NEUMANN:
    # modify Neumann boundary condition on coupling interface, modify weak form correspondingly
    if not boundary_marker:  # there is only 1 Neumann-BC which is at the coupling boundary -> integration over whole boundary
        if coupling_expression.is_scalar_valued():
            # this term has to be added to weak form to add a Neumann BC (see e.g. p.
            # 83ff Langtangen, Hans Petter, and Anders Logg. "Solving PDEs in Python
            # The FEniCS Tutorial Volume I." (2016).)
            F += v * coupling_expression * ds
        elif coupling_expression.is_vector_valued():
            normal = FacetNormal(mesh)
            F += -v * dot(normal, coupling_expression) * ds
        else:
            raise Exception("invalid!")
    else:  # For multiple Neumann BCs integration should only be performed over the respective domain.
        # TODO: fix the problem here
        raise Exception("Boundary markers are not implemented yet")
        # return dot(coupling_bc_expression, v) * dss(boundary_marker)

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")
t = 0

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

# mark mesh w.r.t ranks
mesh_rank = MeshFunction("size_t", mesh, mesh.topology().dim())
if problem is ProblemType.NEUMANN:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 4)
else:
    mesh_rank.set_all(MPI.rank(MPI.comm_world) + 0)
mesh_rank.rename("myRank", "")

# Generating output files
temperature_out = File("output/%s.pvd" % precice.get_participant_name())
ref_out = File("output/ref%s.pvd" % precice.get_participant_name())
error_out = File("output/error%s.pvd" % precice.get_participant_name())
ranks = File("output/ranks%s.pvd" % precice.get_participant_name())

# output solution and reference solution at t=0, n=0
n = 0
print('output u^%d and u_ref^%d' % (n, n))
temperature_out << (u_n, t)
ref_out << u_ref
ranks << mesh_rank

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

flux = Function(V_g)
flux.rename("Flux", "")

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)

    precice_dt = precice.get_max_time_step_size()
    dt.assign(np.min([fenics_dt, precice_dt]))

    # Dirichlet BC and RHS need to point to end of current timestep
    u_D.t = t + float(dt)
    f.t = t + float(dt)

    # Coupling BC needs to point to end of current timestep
    read_data = precice.read_data(dt)

    if problem is ProblemType.DIRICHLET and (domain_part is DomainPart.CIRCULAR or domain_part is DomainPart.RECTANGLE):
        # We have to data for an arbitrary point that is not on the circle, to obtain exact solution.
        # See https://github.com/precice/fenics-adapter/issues/113 for details.
        read_data[(0, 0)] = u_D(0, 0)

    # Update the coupling expression with the new read data
    precice.update_coupling_expression(coupling_expression, read_data)

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    solve(a == L, u_np1, bcs)

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        determine_gradient(V_g, u_np1, flux)
        precice.write_data(flux)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write_data(u_np1)

    precice.advance(dt(0))

    # roll back to checkpoint
    if precice.requires_reading_checkpoint():
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:  # update solution
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        print('n = %d, t = %.2f: L2 error on domain = %.3g' % (n, t, error))
        # output solution and reference solution at t_n+1
        print('output u^%d and u_ref^%d' % (n, n))
        temperature_out << (u_n, t)
        ref_out << u_ref
        error_out << error_pointwise

# Hold plot
precice.finalize()
