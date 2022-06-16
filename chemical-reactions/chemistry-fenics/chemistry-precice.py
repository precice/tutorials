from fenics import *
from mshr import *
import fenicsprecice
import numpy as np
import csv
from mpi4py import MPI


default_dt = 1.0 # time step size

domain = Rectangle(Point(0,0), Point(2.2, 0.41)) - Circle(Point(0.2, 0.2), 0.05)
mesh = generate_mesh(domain, 64)
normal = FacetNormal(mesh)

# Three dimensional vector for three species
P = FiniteElement('P', triangle, 1)
V = FunctionSpace(mesh, MixedElement([P, P, P]))
W = FunctionSpace(mesh, MixedElement([P, P])) # Velocity field is 2D

# No BCs: only homgenuous Von Neumann. It means nothing leaves the domain
# Von neumann BCs give rise to surface loads, which are simply zero here.

class CouplingDomain(SubDomain):
    def inside(self, x, on_boundary):
        return True

# Initialize preCICE
precice = fenicsprecice.Adapter(adapter_config_filename="chemistry-config.json")
precice_dt = precice.initialize(coupling_subdomain=CouplingDomain(), read_function_space=W)

flow_expr = precice.create_coupling_expression()

dt = np.min([default_dt, precice_dt])

# Create used terms (_n suffix means "at previous step")
u_ = Function(V)
u_n = Function(V)
v = TestFunction(V)
u_1, u_2, u_3 = split(u_)
v_1, v_2, v_3 = split(v)
flow = Function(W)
# Formulate the problem

f = Expression(('pow(x[0]-0.1,2)+pow(x[1]-0.1,2)<0.05*0.05 ? 0.1 : 0',
                'pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.1 : 0'), degree=1)

k = Constant(1. / dt)

r = Constant(10.0) # Reaction speed
diff = Constant(0.01) # Diffusivity

# Velocity field when preCICE is not used
# flow = Expression(("vmax * 4*x[1]*(1-x[1])", "0"), degree=2, vmax=0.5)

u_n1, u_n2, u_n3 = split(u_n)
F = k * (u_1 - u_n1) * v_1 * dx + 0.5 * dot(flow, grad(u_1 + u_n1)) * v_1 * dx + diff * dot(grad(u_1), grad(v_1)) * dx - dot(f[0], v_1) * dx + r * u_1 * u_2 * v_1 * dx + \
    k * (u_2 - u_n2) * v_2 * dx + 0.5 * dot(flow, grad(u_2 + u_n2)) * v_2 * dx + diff * dot(grad(u_2), grad(v_2)) * dx - dot(f[1], v_2) * dx + r * u_1 * u_2 * v_2 * dx + \
    k * (u_3 - u_n3) * v_3 * dx + 0.5 * dot(flow, grad(u_3 + u_n3)) * v_3 * dx + diff * dot(grad(u_3), grad(v_3)) * dx                       - r * u_1 * u_2 * v_3 * dx + r * u_3 * v_3 * dx



t = 0
subfolder = 'out_fine_fine_1_order'
vtkfileA = File(subfolder + '/chemical_A.pvd')
vtkfileB = File(subfolder + '/chemical_B.pvd')
vtkfileC = File(subfolder + '/chemical_C.pvd')
vtkfileFlow = File(subfolder + '/chemical_fluid_read.pvd')

# CSV file to keep track of integrals (i.e. total amount of A, B, C)
#with open(subfolder + '/chemical_out.csv', 'w', newline='') as csvfile:

if MPI.COMM_WORLD.rank == 0:
    csvfile = open(subfolder + '/chemical_out.csv', 'w', newline='')
    writer = csv.writer(csvfile, delimiter=' ', quotechar='|',
                      quoting=csv.QUOTE_MINIMAL)

  # No implicit coupling
while precice.is_coupling_ongoing():

    read_data = precice.read_data()
    precice.update_coupling_expression(flow_expr, read_data)
    flow.interpolate(flow_expr)
    # If we add writing, do it here

    dt = np.min([default_dt, precice_dt])
    k.assign(1. / dt)

    t += dt
    print("Time {:.3g}".format(t))
    solve(F == 0, u_)
    u_n.assign(u_)

    u_A, u_B, u_C = u_.split()
    u_A.rename('data', 'data')
    u_B.rename('data', 'data')
    u_C.rename('data', 'data')

    total_A = assemble(u_A * dx)
    total_B = assemble(u_B * dx)
    total_C = assemble(u_C * dx)
    if MPI.COMM_WORLD.rank == 0:
        print(total_A, total_B, total_C)
        writer.writerow([t, total_A, total_B, total_C])

    vtkfileA << u_A, t
    vtkfileB << u_B, t
    vtkfileC << u_C, t
    vtkfileFlow << flow, t

    precice_dt = precice.advance(dt)

precice.finalize()
