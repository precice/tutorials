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
    File, solve, lhs, rhs, grad, inner, dot, dx, ds, interpolate, VectorFunctionSpace, MeshFunction, MPI
from fenicsprecice import Adapter
from ufl_legacy import MixedElement, split

from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
import argparse
import numpy as np
from problem_setup import get_geometry
import dolfin
from dolfin import FacetNormal, dot, project
import sympy as sp
from utils.ButcherTableaux import RadauIIA, BackwardEuler, LobattoIIIC
from utils.high_order_setup import getVariationalProblem, time_derivative
import utils.utils as utl


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


parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
command_group = parser.add_mutually_exclusive_group(required=True)
command_group.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest="dirichlet",
                           action="store_true")
command_group.add_argument("-n", "--neumann", help="create a neumann problem", dest="neumann", action="store_true")
parser.add_argument("-e", "--error-tol", help="set error tolerance", type=float, default=10**-8,)

args = parser.parse_args()

# define problem setup
# get domain part and according problem type
if args.dirichlet and not args.neumann:
    problem = ProblemType.DIRICHLET
    domain_part = DomainPart.LEFT
elif args.neumann and not args.dirichlet:
    problem = ProblemType.NEUMANN
    domain_part = DomainPart.RIGHT

# get computation domain
mesh, coupling_boundary, remaining_boundary = get_geometry(domain_part)

# Define function space using mesh
V = FunctionSpace(mesh, 'P', 2)
V_g = VectorFunctionSpace(mesh, 'P', 1)
W = V_g.sub(0).collapse()

# manufactured solution
alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
# create sympy expression of function to derive the required forms of the equation
x_ ,y_ ,t_ = sp.symbols(['x[0]','x[1]','t'])
u_expr = 1+x_*x_+alpha*y_*y_+beta*t_
u_D = Expression(sp.ccode(u_expr), degree=2, t=0)
# initial condition
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

fenics_dt = .1  # time step size
error_tol = args.error_tol
# define flux function if we are on the dirichlet side of the domain
if problem is ProblemType.DIRICHLET:
    # Define flux in x direction
    flux_expr = Expression(sp.ccode(u_expr.diff(x_)), degree=1, alpha=alpha, t=0)
    flux_fun = interpolate(flux_expr, W)

# time stepping setup
# scheme
tsm = BackwardEuler()
# depending on tsm, we define the trial and test function space
if tsm.num_stages == 1:
    Vbig = V
else:
    # for multi-stage RK methods, we need more dimensional function spaces
    mixed = MixedElement(tsm.num_stages * [V.ufl_element()])
    Vbig = FunctionSpace(V.mesh(), mixed)

# precice setup
precice, precice_dt, initial_data = None, 0.0, None
# Initialize the adapter according to the specific participant
if problem is ProblemType.DIRICHLET:
    precice = Adapter(adapter_config_filename="precice-adapter-config-D.json")
    precice.initialize(coupling_boundary, read_function_space=V, write_object=flux_fun)
    precice_dt = precice.get_max_time_step_size()
elif problem is ProblemType.NEUMANN:
    precice = Adapter(adapter_config_filename="precice-adapter-config-N.json")
    u_D_function = interpolate(u_D, V)
    precice.initialize(coupling_boundary, read_function_space=W, write_object=u_D_function)
    precice_dt = precice.get_max_time_step_size()

dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
# trial and test functions
u = TrialFunction(Vbig)
v = TestFunction(Vbig)
# if dim(Vbig)>1, f needs to be stored in an array with different time stamps,
# because in each stage, of an RK method it is evaluated at a different time
f = tsm.num_stages * [None]
# derive f from sol_expr: f= δu/δt - ∇u
f_expr = u_expr.diff(t_)-u_expr.diff(x_).diff(x_)-u_expr.diff(y_).diff(y_)
for i in range(tsm.num_stages):
    f[i] = Expression(sp.ccode(f_expr), degree=2, t=0)
    f[i].t = tsm.c[i] * float(dt)  # initial time assumed to be 0
# get variational form of the problem
# F = u * v * dx + dt * dot(grad(u), grad(v)) * dx - (u_n + dt * f) * v * dx
F = getVariationalProblem(v=v, initialCondition=u_n, dt=dt, f=f, tsm=tsm, k=u)
a = lhs(F)
L = rhs(F)


# boundary conditions

# we define variational form for each RK stage. According to [Irksome tutorial], All of those stages
# represent time derivatives of the actual solution. Thus, we need time derivatives as boundary conditions

# get time derivative of u
du_dt = time_derivative(u_expr, tsm, dt,t_)
# set up boundary conditions
bc = []
# Create for each dimension of Vbig a coupling expression for either the time derivatives (Dirichlet side)
#  or Neumann side (no time derivatives required as they are enforced with changing F)
coupling_expressions = [precice.create_coupling_expression()] * tsm.num_stages
# for the boundary which is not the coupling boundary, we can just use the boundary conditions as usual
# each stage needs a boundary condition
if tsm.num_stages > 1:
    if problem is ProblemType.DIRICHLET:
        for i in range(tsm.num_stages):
            bc.append(DirichletBC(Vbig.sub(i), du_dt[i], remaining_boundary))
            bc.append(DirichletBC(Vbig.sub(i), coupling_expressions[i], coupling_boundary))
    else:
        vs = split(v)
        for i in range(tsm.num_stages):
            bc.append(DirichletBC(Vbig.sub(i), du_dt[i], remaining_boundary))
            F += vs[i] * coupling_expressions[i] * dolfin.ds
else:
    if problem is ProblemType.DIRICHLET:
        bc.append(DirichletBC(Vbig, du_dt[0], remaining_boundary))
        bc.append(DirichletBC(Vbig, coupling_expressions[0], coupling_boundary))
    else:
        bc.append(DirichletBC(Vbig, du_dt[0], remaining_boundary))
        F += v * coupling_expressions[0] * dolfin.ds


# get lhs and rhs of variational form
a, L = lhs(F), rhs(F)

# we additionally need the values of the stages from the RK method.
# they are stored here
k = Function(Vbig)

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

if problem is ProblemType.DIRICHLET:
    flux = Function(V_g)
    flux.rename("Heat-Flux", "")

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)

    precice_dt = precice.get_max_time_step_size()
    dt.assign(np.min([fenics_dt, precice_dt]))

#    # Dirichlet BC and RHS need to point to end of current timestep
    u_D.t = t + float(dt)
#    f.t = t + float(dt)

    # update boundary conditions
    for i in range(tsm.num_stages):
        f[i].t = t + tsm.c[i] * float(dt)
        du_dt[i].t = t + tsm.c[i] * float(dt)

    # boundary conditions of the coupling boundary needs to be updated as well

    # only dirichlet boundaries need time derivatives
    if problem is ProblemType.DIRICHLET:
        # approximate the function which preCICE uses with BSplines
        bsplns = utl.b_splines(precice, 2, float(dt))
        # get first derivative
        bsplns_der = {}
        for ki in bsplns.keys():
            bsplns_der[ki] = bsplns[ki].derivative(1)

        # preCICE must read num_stages times at respective time for each stage
        for i in range(tsm.num_stages):
            # values of derivative of current time
            val = {}
            for ki in bsplns_der.keys():
                val[ki] = bsplns_der[ki](tsm.c[i]*float(dt))
            precice.update_coupling_expression(coupling_expressions[i], val)
    else:
        # Neumann boundaries just need temperature flux
        for i in range(tsm.num_stages):
            precice.update_coupling_expression(coupling_expressions[i], precice.read_data(tsm.c[i]*dt))

    # getting the solution of the current time step

    # instead of directly solving for u, we look for the values of k.
    # with those we can assemble the solution and thus get the discrete evolution
    solve(a == L, k, bc)

    # now we need to add up the stages k according to the time stepping scheme
    # -> assembly of discrete evolution
    if tsm.num_stages == 1:
        u_np1 = u_n + dt * tsm.b[0] * k
    else:
        u_np1 = u_n
        for i in range(tsm.num_stages):
            u_np1 = u_np1 + dt * tsm.b[i] * k.sub(i)

    # u_sol is in function space V and not Vbig -> project u_np1 to V
    u_sol = project(u_np1, V)

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        determine_gradient(V_g, u_sol, flux)
        flux_x = interpolate(flux.sub(0), W)
        precice.write_data(flux_x)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write_data(u_sol)

    precice.advance(dt)
    precice_dt = precice.get_max_time_step_size()

    # roll back to checkpoint
    if precice.requires_reading_checkpoint():
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:  # update solution
        u_n.assign(u_sol)
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
