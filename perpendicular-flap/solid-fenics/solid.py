# Import required libs
from fenics import Constant, Function, AutoSubDomain, RectangleMesh, VectorFunctionSpace, interpolate, \
    TrialFunction, TestFunction, Point, Expression, DirichletBC, nabla_grad, project, \
    Identity, inner, dx, ds, sym, grad, lhs, rhs, dot, File, solve, PointSource, assemble_system
from ufl import nabla_div
import numpy as np
import matplotlib.pyplot as plt
from fenicsprecice import Adapter
from enum import Enum


# define the two kinds of boundary: clamped and coupling Neumann Boundary
def clamped_boundary(x, on_boundary):
    return on_boundary and abs(x[1]) < tol


def neumann_boundary(x, on_boundary):
    """
    determines whether a node is on the coupling boundary

    """
    return on_boundary and ((abs(x[1] - 1) < tol) or abs(abs(x[0]) - W / 2) < tol)


# Geometry and material properties
dim = 2  # number of dimensions
H = 1
W = 0.1
rho = 3000
E = 4000000
nu = 0.3

mu = Constant(E / (2.0 * (1.0 + nu)))

lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# create Mesh
n_x_Direction = 4
n_y_Direction = 26
mesh = RectangleMesh(Point(-W / 2, 0), Point(W / 2, H), n_x_Direction, n_y_Direction)

h = Constant(H / n_y_Direction)

# create Function Space
V = VectorFunctionSpace(mesh, 'P', 2)

# BCs
tol = 1E-14

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u_np1 = Function(V)
saved_u_old = Function(V)

# function known from previous timestep
u_n = Function(V)
v_n = Function(V)
a_n = Function(V)


f_N_function = interpolate(Expression(("1", "0"), degree=1), V)
u_function = interpolate(Expression(("0", "0"), degree=1), V)

coupling_boundary = AutoSubDomain(neumann_boundary)
fixed_boundary = AutoSubDomain(clamped_boundary)

precice = Adapter(adapter_config_filename="precice-adapter-config-fsi-s.json")

# Initialize the coupling interface
precice_dt = precice.initialize(coupling_boundary, read_function_space=V, write_object=V, fixed_boundary=fixed_boundary)

fenics_dt = precice_dt  # if fenics_dt == precice_dt, no subcycling is applied
# fenics_dt = 0.02  # if fenics_dt < precice_dt, subcycling is applied
dt = Constant(np.min([precice_dt, fenics_dt]))

# clamp the beam at the bottom
bc = DirichletBC(V, Constant((0, 0)), fixed_boundary)

# alpha method parameters
alpha_m = Constant(0)
alpha_f = Constant(0)
gamma = Constant(0.5 + alpha_f - alpha_m)
beta = Constant((gamma + 0.5)**2 / 4.)


# Define strain
def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)


# Define Stress tensor
def sigma(u):
    return lambda_ * nabla_div(u) * Identity(dim) + 2 * mu * epsilon(u)


# Define Mass form
def m(u, v):
    return rho * inner(u, v) * dx


# Elastic stiffness form
def k(u, v):
    return inner(sigma(u), sym(grad(v))) * dx


# External Work
def Wext(u_):
    return dot(u_, p) * ds


# Update functions

# Update acceleration
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)

    return ((u - u_old - dt_ * v_old) / beta / dt_ ** 2
            - (1 - 2 * beta_) / 2 / beta_ * a_old)


# Update velocity
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)

    return v_old + dt_ * ((1 - gamma_) * a_old + gamma_ * a)


def update_fields(u, u_old, v_old, a_old):
    """Update all fields at the end of a timestep."""

    u_vec, u0_vec = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()

    # call update functions
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # assign u->u_old
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()


def avg(x_old, x_new, alpha):
    return alpha * x_old + (1 - alpha) * x_new


# residual
a_np1 = update_a(du, u_n, v_n, a_n, ufl=True)
v_np1 = update_v(a_np1, u_n, v_n, a_n, ufl=True)

res = m(avg(a_n, a_np1, alpha_m), v) + k(avg(u_n, du, alpha_f), v)

a_form = lhs(res)
L_form = rhs(res)

# parameters for Time-Stepping
t = 0.0
n = 0
E_ext = 0

displacement_out = File("Solid/FSI-S/u_fsi.pvd")

u_n.rename("Displacement", "")
u_np1.rename("Displacement", "")
displacement_out << u_n

while precice.is_coupling_ongoing():

    if precice.is_action_required(precice.action_write_iteration_checkpoint()):  # write checkpoint
        precice.store_checkpoint(u_n, t, n)

    # read data from preCICE and get a new coupling expression
    read_data = precice.read_data()

    # Update the point sources on the coupling boundary with the new read data
    Forces_x, Forces_y = precice.get_point_sources(read_data)

    A, b = assemble_system(a_form, L_form, bc)

    b_forces = b.copy()  # b is the same for every iteration, only forces change

    for ps in Forces_x:
        ps.apply(b_forces)
    for ps in Forces_y:
        ps.apply(b_forces)

    assert (b is not b_forces)
    solve(A, u_np1.vector(), b_forces)

    dt = Constant(np.min([precice_dt, fenics_dt]))

    # Write new displacements to preCICE
    precice.write_data(u_np1)

    # Call to advance coupling, also returns the optimum time step value
    precice_dt = precice.advance(dt(0))

    # Either revert to old step if timestep has not converged or move to next timestep
    if precice.is_action_required(precice.action_read_iteration_checkpoint()):  # roll back to checkpoint
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        update_fields(u_np1, saved_u_old, v_n, a_n)
        if n % 10 == 0:
            displacement_out << (u_n, t)

# Plot tip displacement evolution
displacement_out << u_n

precice.finalize()
