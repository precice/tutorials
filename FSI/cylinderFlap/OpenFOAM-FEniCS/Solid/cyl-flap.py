# Import required libs
from fenics import Constant, Function, AutoSubDomain, RectangleMesh, VectorFunctionSpace, interpolate, \
    TrialFunction, TestFunction, Point, Expression, DirichletBC, nabla_grad, \
    Identity, inner, dx, ds, sym, grad, lhs, rhs, dot, File, solve, PointSource, assemble_system
import dolfin

from ufl import nabla_div
import numpy as np
import matplotlib.pyplot as plt
from fenicsadapter import Adapter
from enum import Enum

# Beam geometry
dim = 2  # number of dimensions
L = 0.35  # length
H = 0.02  # height

y_bottom = 0.2 - 0.5 * H  # y coordinate of bottom surface of beam
y_top = y_bottom + H  # y coordinate of top surface of beam
x_left = 0.25  # x coordinate of left surface of beam
x_right = x_left + L  # x coordinate of right surface of beam


# define the two inside functions to determine the boundaries: clamped Dirichlet and coupling Neumann Boundary
def left_boundary(x, on_boundary):
    """
    inside-function for the clamped Dirichlet Boundary.

    Apply Dirichlet boundary on left part of the beam.
    """
    return on_boundary and abs(x[0] - x_left) < tol  # left boundary of beam


def remaining_boundary(x, on_boundary):
    """
    inside-function for the coupling Neumann Boundary.

    Apply Neumann boundary on top, bottom and right (=remaining) part of the beam.
    """
    return on_boundary and (abs(x[1] - y_top) < tol or  # top boundary of beam
                            abs(x[1] - y_bottom) < tol or  # bottom boundary of beam
                            abs(x[0] - x_right) < tol)  # right boundary of beam


# Numerical properties
tol = 1E-14

# Beam material properties
rho = 1000  # density
E = 5600000.0  # Young's modulus
nu = 0.4  # Poisson's ratio
lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))  # first Lame constant
mu = Constant(E / (2.0 * (1.0 + nu)))  # second Lame constant

# create Mesh
n_x_Direction = 20  # DoFs in x direction
n_y_Direction = 4  # DoFs in y direction
mesh = RectangleMesh(Point(x_left, y_bottom), Point(x_right, y_top), n_x_Direction, n_y_Direction)

# create Function Space
V = VectorFunctionSpace(mesh, 'P', 2)

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

# displacement fields
u_np1 = Function(V)
saved_u_old = Function(V)

# function known from previous timestep
u_n = Function(V)
v_n = Function(V)
a_n = Function(V)

# initial value for force and displacement field
f_N_function = interpolate(Expression(("1", "0"), degree=1), V)
u_function = interpolate(Expression(("0", "0"), degree=1), V)

# define coupling boundary
coupling_boundary = AutoSubDomain(remaining_boundary)

adapter_config_filename = "precice-adapter-config-fsi-s.json"

# create Adapter
precice = Adapter(adapter_config_filename)

# create subdomains used by the adapter
clamped_boundary_domain = AutoSubDomain(left_boundary)
force_boundary = AutoSubDomain(remaining_boundary)

# Initialize the coupling interface
precice_dt = precice.initialize(coupling_boundary, mesh, f_N_function, V, dim)

fenics_dt = precice_dt  # if fenics_dt == precice_dt, no subcycling is applied
# fenics_dt = 0.02  # if fenics_dt < precice_dt, subcycling is applied
dt = Constant(np.min([precice_dt, fenics_dt]))

# generalized alpha method (time stepping) parameters
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma = Constant(0.5 + alpha_f - alpha_m)
beta = Constant((gamma + 0.5) ** 2 * 0.25)

# clamp (u == 0) the beam at the left
bc = DirichletBC(V, Constant((0, 0)), left_boundary)


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


def update_acceleration(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)

    return (u - u_old - dt_ * v_old) / beta / dt_ ** 2 - .5 * (1 - 2 * beta_) / beta_ * a_old


def update_velocity(a, u_old, v_old, a_old, ufl=True):
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
    a_vec = update_acceleration(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_velocity(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # assign u->u_old
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()


def avg(x_old, x_new, alpha):
    return alpha * x_old + (1 - alpha) * x_new


# residual
a_np1 = update_acceleration(du, u_n, v_n, a_n, ufl=True)
v_np1 = update_velocity(a_np1, u_n, v_n, a_n, ufl=True)

res = m(avg(a_n, a_np1, alpha_m), v) + k(avg(u_n, du, alpha_f), v)

Forces_x, Forces_y = precice.create_point_sources(clamped_boundary_domain)

a_form = lhs(res)
L_form = rhs(res)

# Prepare for time-stepping
t = 0.0
n = 0
time = []
u_tip = []
time.append(0.0)
u_tip.append(0.0)
E_ext = 0

displacement_out = File("Solid/FSI-S/u_fsi.pvd")

u_n.rename("Displacement", "")
u_np1.rename("Displacement", "")
displacement_out << u_n

# time loop for coupling
while precice.is_coupling_ongoing():

    if precice.is_action_required(precice.action_write_checkpoint()):  # write checkpoint
        precice.store_checkpoint(u_n, t, n)

    # read data from preCICE and get a new coupling expression
    read_data = precice.read_data()

    # Update the point sources on the coupling boundary with the new read data
    Forces_x, Forces_y = precice.update_point_sources(read_data)

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
    if precice.is_action_required(precice.action_read_checkpoint()):  # roll back to checkpoint
        u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_cp)
        t = t_cp
        n = n_cp
    else:
        u_n.assign(u_np1)
        t += dt
        n += 1

    if precice.is_time_window_complete():
        update_fields(u_np1, saved_u_old, v_n, a_n)
        if n % 20 == 0:
            displacement_out << (u_n, t)

        u_tip.append(u_n(0.6, 0.2)[1])
        time.append(t)

# Plot tip displacement evolution
displacement_out << u_n
plt.figure()
plt.plot(time, u_tip)
plt.xlabel("Time")
plt.ylabel("Tip displacement")
plt.show()

precice.finalize()
