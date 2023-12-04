"""
The core functionality of weak_lhs and get_variational_problem were taken from:
    https://github.com/NikoWul/FenicsIrksome
"""

import numpy as np
import scipy.interpolate as interpol
from fenics import *
from sympy import ccode, diff


def b_splines(precice, degree, dt):
    '''
    :param degree: Degree of BSpline
    :param precice: Object of preCICE
    :param dt: current time step size

    Determines an approximation of the waveform from preCICE with BSplines. \n
    The function assumes an equal set of nodes for the whole time interval!
    '''
    # use equidistant samples
    # documentation is not clear about how many samples are required.
    no_samples = 1 + degree
    nodes = np.linspace(0, dt, no_samples)
    weights = np.array([None] * len(nodes))
    b_splns = {}
    key_ref = precice.read_data(0).keys()

    i = 0
    for node in nodes:
        val = precice.read_data(node)
        # check if every key set is equal
        assert(key_ref == val.keys())
        weights[i] = val
        i += 1

    for k in key_ref:
        weights_k = [None] * len(nodes)
        for i in range(no_samples):
            weights_k[i] = weights[i][k]
        b1, b2, b3 = interpol.splrep(nodes, weights_k, s=0, k=degree)
        b_splns[k] = interpol.BSpline(b1, b2, b3, extrapolate=False)
    return b_splns


def weak_lhs(u, v, k):
    """
    Defines the weak left hand side of a given problem
    :param u: trial function
    :param v: test function
    :param k: function of function space for stages of RK methods
    """
    return inner(k, v) * dx + inner(grad(u), grad(v)) * dx


def get_variational_problem(v, k, tsm, f, dt, initial_condition):
    """
    Following https://doi.org/10.1145/3466168, this function creates for each stage of the time stepping scheme
    (Runge-Kutta methods)
    ``tsm``, a weak formulation of a given equation.
    As the stages of Runge-Kutta methods are the time derivatives at different times of the actual solution,
    the solution of the variational problem this function returns, are the time derivatives at the respective stage times.
    It is therefore necessary to assemble after solving the variational form returned by this function
    the discrete evolution according to the time stepping scheme which is used.\n
    **Note that this approach requires the time derivative of Dirichlet boundary conditions.
    Neumann and Robin BCs do not need to be changed.**
    :param v: test function
    :param k: trial function
    :param tsm: time stepping method which should be used
    :param f: rhs of the problem
    :param dt: time step (must not be a float but a fenics expression! e.g. Constant(0))
    :param initial_condition: function which defined the inital value/initial condition
    :return: returns variational problem (F=0) for the given parameters
    """
    num_stages = tsm.num_stages
    ks = split(k)
    vs = split(v)
    u = num_stages * [None]
    for i in range(num_stages):
        uhelp = initial_condition
        for j in range(num_stages):
            uhelp = uhelp + tsm.A[i][j] * dt * ks[j]
        u[i] = uhelp
    rh = 0
    for i in range(num_stages):
        rh = rh + f[i] * vs[i] * dx

    # Assemble weak form from lhs and rhs
    F = 0
    for i in range(num_stages):
        # cf. https://doi.org/10.1145/3466168 p.5 equation 14
        F = F + weak_lhs(v=vs[i], u=u[i], k=ks[i])
    F = F - rh
    return F
