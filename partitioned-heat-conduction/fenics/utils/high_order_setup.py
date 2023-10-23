import sys

from fenics import *
from sympy import ccode, diff


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
        # initial time assumed to be 0
        for j in range(i - 1):
            du_dt[i].t = du_dt[i].t + tsm.c[j] * dt
    return du_dt
