import numpy as np
from numpy.linalg import eig


class SpringLeft:
    k = 4 * np.pi**2


class SpringMiddle:
    k = 16 * (np.pi**2)


class SpringRight:
    k = 4 * np.pi**2


class MassLeft:
    # mass
    m = 1

    # initial conditions (the way how we currently compute the analytical
    # solution allows arbitrary u0, but requires v0 = 0)
    u0 = 1.0
    v0 = 0.0

    u_analytical, v_analytical = None, None  # will be defined below


class MassRight:
    # mass
    m = 1

    # initial conditions (the way how we currently compute the analytical
    # solution allows arbitrary u0, but requires v0 = 0)
    u0 = 0.0
    v0 = 0.0

    u_analytical, v_analytical = None, None  # will be defined below


# Mass matrix
M = np.array([
    [MassLeft.m, 0],
    [0, MassRight.m]
])
# Stiffness matrix
K = np.array([
    [SpringLeft.k + SpringMiddle.k, -SpringMiddle.k],
    [-SpringMiddle.k, SpringRight.k + SpringMiddle.k]
])

# system:
# m ddu + k u = f
# compute analytical solution from eigenvalue ansatz

eigenvalues, eigenvectors = eig(K)
omega = np.sqrt(eigenvalues)
A, B = eigenvectors

c = np.linalg.solve(eigenvectors, [MassLeft.u0, MassRight.u0])

MassLeft.u_analytical = lambda t: c[0] * A[0] * np.cos(omega[0] * t) + c[1] * A[1] * np.cos(omega[1] * t)
MassLeft.v_analytical = lambda t: -c[0] * A[0] * omega[0] * \
    np.sin(omega[0] * t) - c[1] * A[1] * omega[1] * np.sin(omega[1] * t)

MassRight.u_analytical = lambda t: c[0] * B[0] * np.cos(omega[0] * t) + c[1] * B[1] * np.cos(omega[1] * t)
MassRight.v_analytical = lambda t: -c[0] * B[0] * omega[0] * \
    np.sin(omega[0] * t) - c[1] * B[1] * omega[1] * np.sin(omega[1] * t)
