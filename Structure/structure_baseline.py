"""This example is based on implementation of linear elasticity solver from "Langtangen, Hans Petter, and Anders Logg.
Solving PDEs in Python: The FEniCS Tutorial I. Springer International Publishing, 2016. p.50"

It has been updated to the comply with the latest FEniCS version and also modified to 2D scenario (originally 3D).

Program generets files for visulizaiton in ParaView.
"""
from fenics import *
from ufl import nabla_div, nabla_grad
import matplotlib.pyplot as plt

#Scaled varaibles
L = 1 # length of the beam
W = 0.2 # width of square cross section of the beam
mu = 1 # Lame`s elasticity parameters
rho = 1 # density of the beam
delta = W/L
gamma = 0.4*delta**2 # dimensionless variable reflecting ratio of load and shear stress
beta = 1.25 # dimensionless elasticity parameter
lambda_ = beta # Lame`s elasticty parameters
g = gamma # acceleration of gravity

#Mesh and function space
ny = 5 # mesh size in y-direction
nx = 5 # mesh size in x-direction
mesh = RectangleMesh(Point(0, 0), Point(L, W), nx, ny)
V = VectorFunctionSpace(mesh, 'P', 1)

#BC - fixed left edge of the beam
tol = 1E-14
def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol

bc = DirichletBC(V, Constant((0, 0)), clamped_boundary)

#Strain tensor
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    #return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

#Varational problem
u = TrialFunction(V)
d = u.geometric_dimension() # dimension of the space
v = TestFunction(V)
f = Constant((0, -rho*g))
T = Constant((0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

u = Function(V)
solve(a == L, u, bc)

#Plot solution
plt.figure(1)
plot(u, title='Displacement')

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plt.figure(2)
plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
plt.figure(3)
plot(u_magnitude, title='Displacement magnitude')
print('Displacement magnitude u min/max:',
      u_magnitude.vector().min(),
      u_magnitude.vector().max())

#Hold plot
plt.show()

#Save to file
File('output/displacement.pvd') << u
File('output/von_mises.pvd') << von_Mises
File('output/magnitude.pvd') << u_magnitude
