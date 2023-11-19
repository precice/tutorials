import numpy as np
import scipy.interpolate as interpol
from fenics import *
from sympy import ccode, diff


def b_splines(precice, degree, dt):
    '''
    :param degree: Degree of BSpline
    :param precice: Object of precice
    :param dt: current time step size

    Determines an approximation of the waveform from preCICE with BSplines. \n
    The function assumes an equal set of nodes for the whole time interval!
    '''
    # use equidistant samples
    # documentation is not clear about how many samples are required.
    # 2*degree seemingly work fine
    no_samples = 2*degree
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
        b1, b2, b3 = interpol.splrep(nodes, weights_k,s=0, k=degree)
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


def getVariationalProblem(v, k, tsm, f, dt, initialCondition):
    """
    :param v: test function
    :param k: trial function
    :param tsm: time stepping method which should be used
    :param f: f of the problem
    :param dt: time step (must not be a float but a fenics expression! e.g. Constant(0))
    :param initialCondition: function which defined the inital value/initial condition
    :return: returns variational problem (F=0) of the given params from the constructor
    """
    # cf.: Nikolas code in heatConv.py
    num_stages = tsm.num_stages
    ks = split(k)
    vs = split(v)
    u = num_stages * [None]
    for i in range(num_stages):
        uhelp = initialCondition
        for j in range(num_stages):
            uhelp = uhelp + tsm.A[i][j] * dt * ks[j]
        u[i] = uhelp
    rh = 0
    for i in range(num_stages):
        rh = rh + f[i] * vs[i] * dx

    # Assemble weak form from lhs and rhs
    F = 0
    for i in range(num_stages):
        F = F + weak_lhs(v=vs[i], u=u[i], k=ks[i])  # cf. irksome tutorial p.5 eq. 14
    F = F - rh
    return F


def time_derivative(expr, tsm, dt, t):
    """
    :param expr: sympy expression of u
    :param t: sympy symbol for time variable of expr
    :param tsm: time stepping method
    :param dt: length of time step
    :return: expression with the time derivative according to tsm
    """
    # get the time derivative of expr
    du_dt_expr = expr.diff(t)
    du_dt = tsm.num_stages * [None]
    for i in range(tsm.num_stages):
        du_dt[i] = Expression(ccode(du_dt_expr), degree=2, t=0)
        du_dt[i].t = du_dt[i].t + tsm.c[i]*dt(0)
    return du_dt

def determine_gradient(V_g, u, flux):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    :param flux: returns calculated flux into this value
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(grad(u), v) * dx
    solve(a == L, flux)
