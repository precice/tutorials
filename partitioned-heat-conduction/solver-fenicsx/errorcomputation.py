from dolfinx import fem
import numpy as np
from mpi4py import MPI
import ufl


def compute_errors(u_approx, u_ref, total_error_tol=10 ** -9):
    mesh = u_ref.function_space.mesh
    # Compute L2 error and error at nodes
    error_L2 = np.sqrt(mesh.comm.allreduce(fem.assemble_scalar(fem.form((u_approx - u_ref)**2 * ufl.dx)), op=MPI.SUM))
    if mesh.comm.rank == 0:
        print(f"L2-error: {error_L2:.2e}")

    # Compute values at mesh vertices
    error_max = mesh.comm.allreduce(np.max(np.abs(u_approx.x.array - u_ref.x.array)), op=MPI.MAX)
    if mesh.comm.rank == 0:
        print(f"Error_max: {error_max:.2e}")

    assert (error_L2 < total_error_tol)

    return (error_L2, error_max)
