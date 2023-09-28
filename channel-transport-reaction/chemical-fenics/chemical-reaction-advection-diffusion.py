from fenics import Point, FiniteElement, triangle, FunctionSpace, MixedElement, SubDomain, Function, TestFunction, split, Expression, Constant, assemble, solve, dot, grad, dx, File
from mshr import generate_mesh, Rectangle, Circle
import fenicsprecice
import numpy as np
import csv
from mpi4py import MPI
from pathlib import Path

outfolder = Path(__file__).parent / 'output'

default_dt = 1.0  # time step size

domain = Rectangle(Point(0, 0), Point(2.2, 0.41)) - \
    Circle(Point(0.2, 0.2), 0.05)
mesh = generate_mesh(domain, 64)

# Three dimensional vector for three species
P = FiniteElement('P', triangle, 1)
V = FunctionSpace(mesh, MixedElement([P, P, P]))
W = FunctionSpace(mesh, MixedElement([P, P]))  # Velocity field is 2D


class CouplingDomain(SubDomain):
    def inside(self, x, on_boundary):
        return True


# Initialize preCICE
precice = fenicsprecice.Adapter(
    adapter_config_filename=Path(__file__).parent / "chemical-reaction-advection-diffusion.json")
precice.initialize(coupling_subdomain=CouplingDomain(), read_function_space=W)
precice_dt = precice.get_max_time_step_size()

flow_expr = precice.create_coupling_expression()

dt = np.min([default_dt, precice_dt])

# Create used terms (_n suffix means "at previous step")
u_ = Function(V)
u_n = Function(V)
v = TestFunction(V)
u_1, u_2, u_3 = split(u_)
v_1, v_2, v_3 = split(v)
flow = Function(W)
flow_old = Function(W)

# Formulate the problem

f = Expression(
    ('pow(x[0]-0.1,2)+pow(x[1]-0.1,2)<0.05*0.05 ? 0.1 : 0',
     'pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.1 : 0'),
    degree=1)

k = Constant(1. / dt)

r = Constant(10.0)  # Reaction speed
diff = Constant(0.01)  # Diffusivity


average_flow = (flow + flow_old) / 2

u_n1, u_n2, u_n3 = split(u_n)

# Weak form of the PDE. Second order terms (from diffusion) are integrated by parts.
# Integration by parts of laplacian(u) * v * dx yields dot(grad(u), normal) * v * ds - dot(grad(u), grad(v)).
# The surface terms cancels due to the boundary condition dot(grad(u), normal) = 0 (for each u_i).
# Physically speaking, this means that matter doesn't leave the domain
# through diffusion: it only leaves it from the fluid flow.

F = k * (u_1 - u_n1) * v_1 * dx + 0.5 * dot(average_flow,
                                            grad(u_1 + u_n1)) * v_1 * dx + diff * dot(grad(u_1),
                                                                                      grad(v_1)) * dx - dot(f[0],
                                                                                                            v_1) * dx + r * u_1 * u_2 * v_1 * dx + k * (u_2 - u_n2) * v_2 * dx + 0.5 * dot(average_flow,
                                                                                                                                                                                           grad(u_2 + u_n2)) * v_2 * dx + diff * dot(grad(u_2),
                                                                                                                                                                                                                                     grad(v_2)) * dx - dot(f[1],
                                                                                                                                                                                                                                                           v_2) * dx + r * u_1 * u_2 * v_2 * dx + k * (u_3 - u_n3) * v_3 * dx + 0.5 * dot(average_flow,
                                                                                                                                                                                                                                                                                                                                          grad(u_3 + u_n3)) * v_3 * dx + diff * dot(grad(u_3),
                                                                                                                                                                                                                                                                                                                                                                                    grad(v_3)) * dx - r * u_1 * u_2 * v_3 * dx + r * u_3 * v_3 * dx


t = 0
vtkfileA = File(str(outfolder / 'chemical_A.pvd'))
vtkfileB = File(str(outfolder / 'chemical_B.pvd'))
vtkfileC = File(str(outfolder / 'chemical_C.pvd'))
vtkfileFlow = File(str(outfolder / 'chemical_fluid_read.pvd'))

# CSV file to keep track of integrals (i.e. total amount of A, B, C)
# with open(outfolder / 'chemical_out.csv', 'w', newline='') as csvfile:

if MPI.COMM_WORLD.rank == 0:
    csvfile = open(outfolder / 'chemical_out.csv', 'w', newline='')
    writer = csv.writer(csvfile, delimiter=' ', quotechar='|',
                        quoting=csv.QUOTE_MINIMAL)

  # No implicit coupling
while precice.is_coupling_ongoing():

    precice_dt = precice.get_max_time_step_size()
    dt = np.min([default_dt, precice_dt])
    read_data = precice.read_data(dt)
    precice.update_coupling_expression(flow_expr, read_data)
    flow_old.assign(flow)
    flow.interpolate(flow_expr)
    # If we add writing, do it here

    k.assign(1. / dt)

    t += dt
    print("Time {:.3g}".format(t))
    solve(F == 0, u_)
    u_n.assign(u_)

    u_A, u_B, u_C = u_.split()
    u_A.rename('data', 'data')
    u_B.rename('data', 'data')
    u_C.rename('data', 'data')

    # Compute the integrals. (All ranks must do it to synchronize)
    total_A = assemble(u_A * dx)
    total_B = assemble(u_B * dx)
    total_C = assemble(u_C * dx)
    # Only the main rank writes to the file.
    if MPI.COMM_WORLD.rank == 0:
        print(total_A, total_B, total_C)
        writer.writerow([t, total_A, total_B, total_C])

    vtkfileA << u_A, t
    vtkfileB << u_B, t
    vtkfileC << u_C, t
    vtkfileFlow << flow, t

    precice.advance(dt)

precice.finalize()
