from ufl import dx
from dolfinx.fem import assemble_scalar, form
import numpy as np
from mpi4py import MPI


def compute_errors(u_approx, u_ref, total_error_tol=10 ** -4):
    mesh = u_ref.function_space.mesh

    # compute total L2 error between reference and calculated solution
    error_pointwise = form(((u_approx - u_ref) / u_ref) ** 2 * dx)
    error_total = np.sqrt(mesh.comm.allreduce(assemble_scalar(error_pointwise), MPI.SUM))

    assert (error_total < total_error_tol)

    return error_total
