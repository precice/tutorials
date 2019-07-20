import fenics
from fenics import MPI, Point, FunctionSpace, Expression, interpolate, RectangleMesh

print("Hello World from rank {rank} of {size}.".format(rank=MPI.rank(MPI.comm_world), size=MPI.size(MPI.comm_world)))

nx = 10
ny = 10
alpha = 3  # parameter alpha
beta = 1.3  # parameter beta
y_bottom, y_top = 0, 1
x_left, x_right = 0, 2
x_coupling = 1.5  # x coordinate of coupling interface

p0 = Point(x_left, y_bottom)
p1 = Point(x_coupling, y_top)

mesh = RectangleMesh(p0, p1, nx, ny)
V = FunctionSpace(mesh, 'P', 1)

u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0)
u_D_function = interpolate(u_D, V)

for v in fenics.vertices(mesh):
    print("{rank} of {size}: u_D = {data} @ ({x},{y})".format(rank=MPI.rank(MPI.comm_world),
                                                              size=MPI.size(MPI.comm_world),
                                                              data=u_D_function(v.x(0), v.x(1)),
                                                              x=v.x(0),
                                                              y=v.x(1)))