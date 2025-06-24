"""
This tutorial is an adaption to preCICE of a case introduced in the FEniCS tutorial book:
"Solving PDEs in Python" by Hans Petter Langtangen and Anders Logg.
See https://fenicsproject.org/tutorial/

Instead of loading data to solve chemical reactions, we send the velocity field to preCICE.
"""


from fenics import *
from mshr import *
import fenicsprecice
import numpy as np
from pathlib import Path

outfolder = Path(__file__).parent / 'output'

# time step size (preCICE handles the main loop; if this is smaller than
# preCICE dt, we subcycle)
default_dt = 0.001
mu = 0.001  # Dynamic viscosity
rho = 1  # Density

domain = Rectangle(Point(0, 0), Point(2.2, 0.41)) - \
    Circle(Point(0.2, 0.2), 0.05)
mesh = generate_mesh(domain, 50)
normal = FacetNormal(mesh)

# Expressions for evaluating BC
inflow = "near(x[0], 0)"
outflow = "near(x[0], 2.2)"
walls = "near(x[1], 0) || near(x[1], 0.41)"
cylinder = "on_boundary && x[0] > 0.1 && x[0] < 0.3 && x[1] > 0.1 && x[1] < 0.3"

# Second order elements for velocity and first order for pressure
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

# Define parabolic inflow profile
inflow_profile = ("4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)", "0")

# Boundary conditions
# Imposed velocity in the outlet and free pressure at the outlet. Walls
# are not slipping
bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), inflow)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]


# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_ = Function(V)
p_n = Function(Q)
p_ = Function(Q)

# Define expressions used in variational forms
U = 0.5 * (u_n + u)
n = FacetNormal(mesh)
f = Constant((0, 0))
k = Constant(default_dt)
mu = Constant(mu)
rho = Constant(rho)

# Define symmetric gradient


def epsilon(u):
    return sym(nabla_grad(u))
# Define stress tensor


def sigma(u, p):
    return 2 * mu * epsilon(u) - p * Identity(len(u))


# Define variational problem for step 1
F1 = rho * dot((u - u_n) / k, v) * dx \
    + rho * dot(dot(u_n, nabla_grad(u_n)), v) * dx \
    + inner(sigma(U, p_n), epsilon(v)) * dx \
    + dot(p_n * n, v) * ds - dot(mu * nabla_grad(U) * n, v) * ds \
    - dot(f, v) * dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q)) * dx
L2 = dot(nabla_grad(p_n), nabla_grad(q)) * dx - (1 / k) * div(u_) * q * dx

# Define variational problem for step 3
a3 = dot(u, v) * dx
L3 = dot(u_, v) * dx - k * dot(nabla_grad(p_ - p_n), v) * dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Initialize preCICE

# Coupling subdomain is the entire domain (volumetric coupling)


class CouplingDomain(SubDomain):
    def inside(self, x, on_boundary):
        return True


precice = fenicsprecice.Adapter(adapter_config_filename=Path(__file__).parent / "fluid-config.json")
precice.initialize(coupling_subdomain=CouplingDomain(), write_object=u_)
precice_dt = precice.get_max_time_step_size()

dt = np.min([default_dt, precice_dt])
k.assign(dt)

# No implicit coupling supported, as this is uni-directional coupling
# If needed, implement checkpointing
t = 0
vtkfile = File(str(outfolder / 'chemical_fluid_write.pvd'))

while precice.is_coupling_ongoing():

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')
    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')
    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg', 'sor')

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

    precice.write_data(u_)

    t += dt

    vtkfile << u_

    precice.advance(dt)
    precice_dt = precice.get_max_time_step_size()
    dt = np.min([default_dt, precice_dt])
    k.assign(dt)

precice.finalize()
