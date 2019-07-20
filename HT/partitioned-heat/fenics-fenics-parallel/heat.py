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
from fenics import Function, SubDomain, RectangleMesh, FunctionSpace, Point, Expression, Constant, DirichletBC, \
    TrialFunction, TestFunction, File, solve, plot, lhs, rhs, grad, inner, dot, dx, ds, assemble, interpolate, project, near, MPI
from enum import Enum
import fenicsadapter.core
from errorcomputation import compute_errors
import argparse
import numpy as np


class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem

class Subcyling(Enum):
    """
    Enum defines which kind of subcycling is used
    """
    NONE = 0  # no subcycling, precice_dt == fenics_dt
    MATCHING = 1  # subcycling, where fenics_dt fits into precice_dt, mod(precice_dt, fenics_dt) == 0
    NONMATCHING = 2  # subcycling, where fenics_dt does not fit into precice_dt, mod(precice_dt, fenics_dt) != 0

    # note: the modulo expressions above should be understood in an exact way (no floating point round off problems. For
    # details, see https://stackoverflow.com/questions/14763722/python-modulo-on-floats)


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
        if on_boundary and near(x[0], x_coupling, tol) and not near(x[1], y_bottom, tol) and not near(x[1], y_top, tol):
            return True
        else:
            return False


def fluxes_from_temperature_full_domain(F, V):
    """
    compute flux from weak form (see p.3 in Toselli, Andrea, and Olof Widlund. Domain decomposition methods-algorithms and theory. Vol. 34. Springer Science & Business Media, 2006.)
    :param F: weak form with known u^{n+1}
    :param V: function space
    :param hy: spatial resolution perpendicular to flux direction
    :return:
    """
    print("rank {rank} of {size}: enters fluxes_from_temperature_full_domain".format(rank=MPI.rank(MPI.comm_world), size=MPI.size(MPI.comm_world)))
    fluxes_vector = assemble(F)  # assemble weak form -> evaluate integral
    print("rank {rank} of {size}: first assembly done".format(rank=MPI.rank(MPI.comm_world),
                                                                                     size=MPI.size(MPI.comm_world)))
    v = TestFunction(V)
    myrank = MPI.rank(MPI.comm_world)
    fluxes = MPI.size(MPI.comm_world) * [None]
    fluxes[myrank] = Function(V)  # create function for flux
    area = assemble(v * ds).get_local()
    """
    print("rank {rank} of {size}: second assembly done".format(rank=MPI.rank(MPI.comm_world),
                                                                                     size=MPI.size(MPI.comm_world)))
    print("rank {rank} has to visit {n_vertices} vertices.".format(rank=MPI.rank(MPI.comm_world), n_vertices=fluxes_vector.__len__()))
    for i in range(area.shape[0]):
        print("rank {rank}: fluxes_vector[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i, value=fluxes_vector[i]))
        print("rank {rank}: area[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i,
                                                               value=area[i]))
        print("rank {rank}: fluxes.vector()[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i,
                                                                 value=fluxes[myrank].vector()[i]))
    for i in range(area.shape[0]):
        print("rank {rank} of {size}: visits DoF {i}".format(rank=MPI.rank(MPI.comm_world),
                                                                   size=MPI.size(MPI.comm_world),
                                                             i=i))
        if area[i] != 0:  # put weight from assemble on function
            print("rank {rank} of {size}: case 1".format(rank=MPI.rank(MPI.comm_world),
                                                                                     size=MPI.size(MPI.comm_world)))
            print("rank {rank}: fluxes_vector[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i,
                                                                   value=fluxes_vector[i]))
            print("rank {rank}: area[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i,
                                                          value=area[i]))
            print("rank {rank}: fluxes.vector()[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i,
                                                                     value=fluxes[myrank].vector()[i]))
            fluxes[myrank].vector()[i] = fluxes_vector[i] / area[i]  # scale by surface area
            print("rank {rank}: fluxes.vector()[{i}]={value}".format(rank=MPI.rank(MPI.comm_world), i=i,
                                                                     value=fluxes[myrank].vector()[i]))
            print("rank {rank} of {size}: computation done".format(rank=MPI.rank(MPI.comm_world),
                                                                                     size=MPI.size(MPI.comm_world)))
        else:
            print("rank {rank} of {size}: case 2".format(rank=MPI.rank(MPI.comm_world),
                                                                                     size=MPI.size(MPI.comm_world)))
            assert(abs(fluxes_vector[i]) < 10**-10)  # for non surface parts, we expect zero flux
            fluxes[myrank].vector()[i] = fluxes_vector[i]
        print("rank {rank} of {size}: completed DoF {i}".format(rank=MPI.rank(MPI.comm_world),
                                                             size=MPI.size(MPI.comm_world),
                                                             i=i))
    print("rank {rank} of {size}: waits...".format(rank=MPI.rank(MPI.comm_world),
                                                   size=MPI.size(MPI.comm_world)))
    """
    MPI.barrier(MPI.comm_world)
    print("rank {rank} of {size}: exits fluxes_from_temperature_full_domain".format(rank=MPI.rank(MPI.comm_world),
                                                                                     size=MPI.size(MPI.comm_world)))
    return fluxes_vector#[myrank]


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dirichlet", help="create a dirichlet problem", dest='dirichlet', action='store_true')
parser.add_argument("-n", "--neumann", help="create a neumann problem", dest='neumann', action='store_true')

args = parser.parse_args()


config_file_name = "precice-config.xml"  # TODO should be moved into config, see https://github.com/precice/fenics-adapter/issues/5 , in case file doesnt not exsist open will fail

# coupling parameters
if args.dirichlet:
    problem = ProblemType.DIRICHLET
    solver_name = "HeatDirichlet"
    mesh_name = "DirichletNodes"
    write_data_name = "Flux"
    read_data_name = "Temperature"
if args.neumann:
    problem = ProblemType.NEUMANN
    solver_name = "HeatNeumann"
    mesh_name = "NeumannNodes"
    write_data_name = "Temperature"
    read_data_name = "Flux"
if args.dirichlet and args.neumann:
    raise Exception("you can only choose either a dirichlet problem (option -d) or a neumann problem (option -n)")
if not (args.dirichlet or args.neumann):
    raise Exception("you have to choose either a dirichlet problem (option -d) or a neumann problem (option -n)")

# Create mesh and define function space

nx = 5
ny = 10
subcycle = Subcyling.NONE

assert(MPI.initialized())

if problem is ProblemType.DIRICHLET:
    nx = nx*3

elif problem is ProblemType.NEUMANN:
    ny = 20

# for all scenarios, we assume precice_dt == .1
if subcycle is Subcyling.NONE:
    fenics_dt = .1  # time step size
    error_tol = 10 ** -4  # error low, if we do not subcycle. In theory we would obtain the analytical solution.
    # TODO For reasons, why we currently still have a relatively high error, see milestone https://github.com/precice/fenics-adapter/milestone/1
elif subcycle is Subcyling.MATCHING:
    fenics_dt = .01  # time step size
    error_tol = 10 ** -2  # error increases. If we use subcycling, we cannot assume that we still get the exact solution.
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.
elif subcycle is Subcyling.NONMATCHING:
    fenics_dt = .03  # time step size
    error_tol = 10 ** -1  # error increases. If we use subcycling, we cannot assume that we still get the exact solution.
    # TODO Using waveform relaxation, we should be able to obtain the exact solution here, as well.

alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface

if problem is ProblemType.DIRICHLET:
    p0 = Point(x_left, y_bottom)
    p1 = Point(x_coupling, y_top)
elif problem is ProblemType.NEUMANN:
    p0 = Point(x_coupling, y_bottom)
    p1 = Point(x_right, y_top)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0)
u_D_function = interpolate(u_D, V)
# Define flux in x direction on coupling interface (grad(u_D) in normal direction)
f_N = Expression('2 * x[0]', degree=1)
f_N_function = interpolate(f_N, V)

coupling_boundary = CouplingBoundary()
remaining_boundary = ComplementaryBoundary(coupling_boundary)

bcs = [DirichletBC(V, u_D, remaining_boundary)]
# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")

if problem is ProblemType.DIRICHLET:
    write_field = f_N_function
    read_field = u_D_function
elif problem is ProblemType.NEUMANN:
    write_field = u_D_function
    read_field = f_N_function

print("Hello from rank {rank} of {size}.".format(rank=MPI.rank(MPI.comm_world), size=MPI.size(MPI.comm_world)))

precice = fenicsadapter.core.Adapter(solver_name, MPI.rank(MPI.comm_world), MPI.size(MPI.comm_world))
precice.configure(config_file_name)
precice.set_coupling_mesh(mesh_name, mesh, coupling_boundary)
precice_dt = precice.initialize()

if precice.is_action_required(fenicsadapter.core.action_write_initial_data()):
    precice.write_block_scalar_data(write_data_name, mesh_name, write_field)
    precice.fulfilled_action(fenicsadapter.core.action_write_initial_data())

precice.initialize_data()

if precice.is_read_data_available():
    coupling_expression = precice.read_block_scalar_data(read_data_name, mesh_name)
else:
    coupling_expression = precice.create_coupling_boundary_condition(read_field)

dt = Constant(0)
dt.assign(np.min([fenics_dt, precice_dt]))

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2 * alpha)
F = u * v / dt * dx + dot(grad(u), grad(v)) * dx - (u_n / dt + f) * v * dx

if problem is ProblemType.DIRICHLET:
    # apply Dirichlet boundary condition on coupling interface
    bcs.append(DirichletBC(V, coupling_expression, coupling_boundary))
if problem is ProblemType.NEUMANN:
    # apply Neumann boundary condition on coupling interface, modify weak form correspondingly
    F += coupling_expression * v * ds

a, L = lhs(F), rhs(F)

# Time-stepping
u_np1 = Function(V)
F_known_u = u_np1 * v / dt * dx + dot(grad(u_np1), grad(v)) * dx - (u_n / dt + f) * v * dx
u_np1.rename("Temperature", "")

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

temperature_out = File("out/%s.pvd" % solver_name)
ref_out = File("out/ref%s.pvd" % solver_name)
error_out = File("out/error%s.pvd" % solver_name)

# output solution and reference solution at t=0, n=0
t = 0
n = 0
print('{rank}:output u^{iteration} and u_ref^{iteration}'.format(rank=MPI.rank(MPI.comm_world), iteration=n))
temperature_out << u_n
ref_out << u_ref

error_total, error_pointwise = compute_errors(u_n, u_ref, V)
error_out << error_pointwise

# set t_1 = t_0 + dt, this gives u_D^1
u_D.t = t + dt(0)  # call dt(0) to evaluate FEniCS Constant. Todo: is there a better way?

# write checkpoint
if precice.is_action_required(fenicsadapter.core.action_write_iteration_checkpoint()):
    u_cp = u_n.copy(deepcopy=True)
    t_cp = t
    n_cp = n
    precice.fulfilled_action(fenicsadapter.core.action_write_iteration_checkpoint())

while precice.is_coupling_ongoing():

    # Compute solution u^n+1, use bcs u_D^n+1, u^n and coupling bcs
    print('{rank} of {size}:starts solving'.format(rank=MPI.rank(MPI.comm_world), size=MPI.size(MPI.comm_world)))
    solve(a == L, u_np1, bcs)
    print('{rank}:done solving'.format(rank=MPI.rank(MPI.comm_world)))

    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
        write_field = fluxes_from_temperature_full_domain(F_known_u, V)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem obtains sends temperature on boundary to Dirichlet problem
        write_field = u_np1
    MPI.barrier(MPI.comm_world)
    print('{rank} of {size}:before write'.format(rank=MPI.rank(MPI.comm_world), size=MPI.size(MPI.comm_world)))
    precice.write_block_scalar_data(write_data_name, mesh_name, write_field)
    precice_dt = precice.advance(dt)
    dt.assign(np.min([fenics_dt, precice_dt]))
    read_expression = precice.read_block_scalar_data(read_data_name, mesh_name)
    coupling_expression.update(read_expression)

    # checkpointing
    if precice.is_action_required(fenicsadapter.core.action_read_iteration_checkpoint()):
        # continue FEniCS computation from checkpoint
        u_n.assign(u_cp)  # set u_n to value of checkpoint
        t = t_cp
        n = n_cp
        precice.fulfilled_action(fenicsadapter.core.action_read_iteration_checkpoint())
    else:
        u_n.assign(u_np1)
        t = t + dt
        n = n + 1

    if precice.is_action_required(fenicsadapter.core.action_write_iteration_checkpoint()):
        # continue FEniCS computation with u_np1
        # update checkpoint
        u_cp.assign(u_np1)
        t_cp = t
        n_cp = n
        precice.fulfilled_action(fenicsadapter.core.action_write_iteration_checkpoint())
        # write results and chec
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        print('n = %d, t = %.2f: L2 error on domain = %.3g' % (n, t, error))
        # output solution and reference solution at t_n+1
        print('output u^%d and u_ref^%d' % (n, n))
        temperature_out << u_n
        ref_out << u_ref
        error_out << error_pointwise

    # Update dirichlet BC
    u_D.t = t + dt(0)

# Hold plot
precice.finalize()

assert(MPI.finalized())