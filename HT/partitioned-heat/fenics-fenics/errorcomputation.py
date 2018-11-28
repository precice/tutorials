from fenics import inner, assemble, dx, project

def compute_errors(u_approx, u_ref, V, total_error_tol = 10**-4):
    error = inner(u_ref - u_approx, u_ref - u_approx)/(u_ref * u_ref)  # compute pointwise L2 error
    error_pointwise = project(inner(u_ref - u_approx, u_ref - u_approx)/(u_ref * u_ref), V)  # project onto function space
    error_total = assemble(error_pointwise * dx)  # determine L2 norm to estimate total error
    print(error_pointwise)
    error_pointwise.rename("error", " ")

    assert (error_total < total_error_tol)

    return error_total, error_pointwise
