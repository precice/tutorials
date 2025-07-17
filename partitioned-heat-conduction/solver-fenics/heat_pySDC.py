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


This variant of this tutorial example uses the open-source library pySDC for time-stepping.
pySDC can be installed via `pip install pySDC`.
If you want to use the developer version, the source code repository can be cloned from "https://github.com/Parallel-in-Time/pySDC".
For more information on pySDC, see also "https://parallel-in-time.org/pySDC/"

For information on the partitioned heat conduction problem using pySDC methods see
* "Eppacher, Tobias. Parallel-in-Time Integration with preCICE. Bachelor's thesis at Technical University of Munich, 2024. URL: https://mediatum.ub.tum.de/1755012"

"""

from __future__ import print_function, division
from fenics import Function, FunctionSpace, Expression, Constant, TrialFunction, TestFunction, \
    File, solve, grad, inner, dx, interpolate, VectorFunctionSpace, MeshFunction, MPI
from fenicsprecice import Adapter
from errorcomputation import compute_errors
from my_enums import ProblemType, DomainPart
import argparse
import numpy as np
from problem_setup import get_geometry
import sympy as sp

from pySDC.implementations.sweeper_classes.imex_1st_order_mass import imex_1st_order_mass
from pySDC.implementations.controller_classes.controller_nonMPI import controller_nonMPI
from heat_pySDC_problem_class import fenics_heat_2d


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


def setup_problem(
        function_space,
        coupling_boundary,
        remaining_boundary,
        u_D, forcing_expr,
        coupling_expression,
        precice,
        dt,
        logger_level=30,
        quad_type='LOBATTO',
        num_nodes=4,
        restol=1e-11,
        maxiter=40):

    # Create docstring for this function
    """
    Setup the problem and controller for the heat equation problem.

    Args:
        function_space: FEniCS function space object
        coupling_boundary: FEniCS SubDomain object for the coupling boundary
        remaining_boundary: FEniCS SubDomain object for the remaining boundary
        u_D: FEniCS expression for the manufactured solution
        forcing_expr: FEniCS expression for the forcing term
        coupling_expression: FEniCS expression for the coupling boundary condition
        precice: preCICE-FEniCS adapter object reference
        dt: time step size
        logger_level: logging level
        quad_type: quadrature type
        num_nodes: number of nodes
        restol: residual tolerance
        maxiter: maximum number of iterations

    Returns:
        controller: pySDC controller object
        P: problem object
    """

    # initialize controller parameters
    controller_params = {
        'logger_level': logger_level
    }

    # fill description dictionary for easy instantiation
    description = {
        'problem_class': fenics_heat_2d,
        'problem_params': {
            'function_space': function_space,
            'coupling_boundary': coupling_boundary,
            'remaining_boundary': remaining_boundary,
            'solution_expr': u_D,
            'forcing_term_expr': forcing_expr,
            'precice_ref': precice,
            'coupling_expr': coupling_expression
        },
        'sweeper_class': imex_1st_order_mass,
        'sweeper_params': {
            'quad_type': quad_type,
            'num_nodes': num_nodes,
        },
        'level_params': {
            'restol': restol,
            'dt': dt
        },
        'step_params': {
            'maxiter': maxiter,
        }
    }

    # Controller for time stepping
    controller = controller_nonMPI(num_procs=1, controller_params=controller_params, description=description)

    # Reference to problem class for easy access to exact solution
    P = controller.MS[0].levels[0].prob
    return controller, P


parser = argparse.ArgumentParser(description="Solving heat equation for simple or complex interface case")
parser.add_argument("participantName", help="Name of the solver.", type=str, choices=[p.value for p in ProblemType])
parser.add_argument("-e", "--error-tol", help="set error tolerance", type=float, default=10**-8,)

args = parser.parse_args()
participant_name = args.participantName


# preCICE error tolerance
# Error is bounded by coupling accuracy. In theory we would obtain the analytical solution.
error_tol = args.error_tol

if participant_name == ProblemType.DIRICHLET.value:
    problem = ProblemType.DIRICHLET
    domain_part = DomainPart.LEFT
elif participant_name == ProblemType.NEUMANN.value:
    problem = ProblemType.NEUMANN
    domain_part = DomainPart.RIGHT
    exit("Neumann problem not yet supported with pySDC")

domain_mesh, coupling_boundary, remaining_boundary = get_geometry(domain_part)

# Define function space using mesh
V = FunctionSpace(domain_mesh, 'P', 2)
V_g = VectorFunctionSpace(domain_mesh, 'P', 1)
W = V_g.sub(0).collapse()


# Time step size
# Should be integer fraction of the used time window size
pySDC_dt = 0.1

# Manufactured solution parameters
alpha = 3  # parameter alpha
beta = 1.2  # parameter beta


# create sympy expression of manufactured solution
x_sp, y_sp, t_sp = sp.symbols(['x[0]', 'x[1]', 't'])

temporal_deg = 1  # temporal degree of the manufactured solution
u_D_sp = 1 + x_sp * x_sp + alpha * y_sp * y_sp + beta * (t_sp ** temporal_deg)
u_D = Expression(sp.ccode(u_D_sp), degree=2, alpha=alpha, beta=beta, temporal_deg=temporal_deg, t=0)

u_D_function = interpolate(u_D, V)
f_sp = u_D_sp.diff(t_sp) - u_D_sp.diff(x_sp).diff(x_sp) - u_D_sp.diff(y_sp).diff(y_sp)
forcing_expr = Expression(sp.ccode(f_sp), degree=2, alpha=alpha, beta=beta, t=0)


if problem is ProblemType.DIRICHLET:
    # Define flux in x direction
    f_N = Expression(sp.ccode(u_D_sp.diff(x_sp)), degree=1, alpha=alpha, t=0)
    f_N_function = interpolate(f_N, W)

# Define initial value
u_n = interpolate(u_D, V)
u_n.rename("Temperature", "")


precice, precice_dt, initial_data = None, 0.0, None

# Initialize the adapter according to the specific participant
precice = Adapter(adapter_config_filename="precice-adapter-config.json")

if problem is ProblemType.DIRICHLET:
    precice.initialize(coupling_boundary, read_function_space=V, write_object=f_N_function)
elif problem is ProblemType.NEUMANN:
    precice.initialize(coupling_boundary, read_function_space=W, write_object=u_D_function)

precice_dt = precice.get_max_time_step_size()
dt = Constant(0)
dt.assign(np.min([pySDC_dt, precice_dt]))

coupling_expression = precice.create_coupling_expression()


controller, P = setup_problem(V,
                              coupling_boundary,
                              remaining_boundary,
                              u_D, forcing_expr,
                              coupling_expression,
                              precice,
                              pySDC_dt)


# Time-stepping
u_np1 = Function(V)
u_np1.rename("Temperature", "")
t = 0

# reference solution at t=0
u_ref = interpolate(u_D, V)
u_ref.rename("reference", " ")

# mark mesh w.r.t ranks
mesh_rank = MeshFunction("size_t", domain_mesh, domain_mesh.topology().dim())
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
print("output u^%d and u_ref^%d" % (n, n))
ranks << mesh_rank

error_total, error_pointwise = compute_errors(u_n, u_ref, V)

# create buffer for output. We need this buffer, because we only want to
# write the converged output at the end of the window, but we also want to
# write the samples that are resulting from substeps inside the window
u_write = []
ref_write = []
error_write = []
# copy data to buffer and rename
uu = u_n.copy()
uu.rename("u", "")
u_write.append((uu, t))
uu_ref = u_ref.copy()
uu_ref.rename("u_ref", "")
ref_write.append(uu_ref)
err = error_pointwise.copy()
err.rename("err", "")
error_write.append(err)

if problem is ProblemType.DIRICHLET:
    flux = Function(V_g)
    flux.rename("Heat-Flux", "")

while precice.is_coupling_ongoing():

    # write checkpoint
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)

        # output solution and reference solution at t_n+1 and substeps (read from buffer)
        print('output u^%d and u_ref^%d' % (n, n))
        for sample in u_write:
            temperature_out << sample

        for sample in ref_write:
            ref_out << sample

        for sample in error_write:
            error_out << error_pointwise

    precice_dt = precice.get_max_time_step_size()
    dt.assign(np.min([pySDC_dt, precice_dt]))

    # Update start time of the current preCICE time step within
    # the problem class.
    # This is necessary, that the coupling expression update in the problem class
    # is executed correctly.
    P.t_start = t

    # Retrieve the result at the end of the timestep from pySDC
    uend, _ = controller.run(u_n, t0=t, Tend=t + float(dt))

    # Update the buffer with the new solution values
    u_np1 = uend.values

    # Write data to preCICE according to which problem is being solved
    if problem is ProblemType.DIRICHLET:
        # Dirichlet problem reads temperature and writes flux on boundary to Neumann problem
        determine_gradient(V_g, u_np1, flux)
        flux_x = interpolate(flux.sub(0), W)
        precice.write_data(flux_x)
    elif problem is ProblemType.NEUMANN:
        # Neumann problem reads flux and writes temperature on boundary to Dirichlet problem
        precice.write_data(u_np1)

    precice.advance(dt)
    precice_dt = precice.get_max_time_step_size()

    # roll back to checkpoint
    if precice.requires_reading_checkpoint():
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
        # empty buffer if window has not converged
        u_write = []
        ref_write = []
        error_write = []
    else:  # update solution
        u_n.assign(u_np1)
        t += float(dt)
        n += 1
        # copy data to buffer and rename
        uu = u_n.copy()
        uu.rename("u", "")
        u_write.append((uu, t))
        uu_ref = u_ref.copy()
        uu_ref.rename("u_ref", "")
        ref_write.append(uu_ref)
        err = error_pointwise.copy()
        err.rename("err", "")
        error_write.append(err)

    if precice.is_time_window_complete():
        u_ref = interpolate(u_D, V)
        u_ref.rename("reference", " ")
        error, error_pointwise = compute_errors(u_n, u_ref, V, total_error_tol=error_tol)
        print("n = %d, t = %.2f: L2 error on domain = %.3g" % (n, t, error))


# output solution and reference solution at t_n+1 and substeps (read from buffer)
print("output u^%d and u_ref^%d" % (n, n))
for sample in u_write:
    temperature_out << sample

for sample in ref_write:
    ref_out << sample

for sample in error_write:
    error_out << error_pointwise

# Hold plot
precice.finalize()
