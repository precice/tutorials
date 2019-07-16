from fenics import MPI

print("Hello World from rank {rank} of {size}.".format(rank=MPI.rank(MPI.comm_world), size=MPI.size(MPI.comm_world)))
