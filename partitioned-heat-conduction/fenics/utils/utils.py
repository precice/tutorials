import argparse
import sys
import numpy as np

from fenics import dx, TrialFunction, TestFunction, solve, inner, grad
import scipy.interpolate as interpol

def b_splines(precice, degree, dt):
    '''
    :param degree: Degree of BSpline
    :param precice: Object of precice
    :param dt: current time step size
    '''
    # use equidistant samples
    # apparently you need 2k+2 samples for degree k
    x_dist = float(dt) / (2 * degree + 2)
    nodes = []
    weights = []
    b_splines = {}
    for i in range(2 * degree + 3):
        if i==0:
            nodes.append(0)
            weights.append(precice.read_data(0))
        elif i==2*degree+2:
            nodes.append(dt)
            weights.append(precice.read_data(dt))
        else:
            nodes.append(x_dist * i)
            weights.append(precice.read_data(x_dist * i))
    for k in weights[0].keys():
        weights_k = []
        for i in range(2 * degree + 3):
            weights_k.append(weights[i][k])
        b1,b2,b3 = interpol.splrep(nodes, weights_k,s=0, k=degree)
        b_splines[k] = interpol.BSpline(b1, b2, b3, extrapolate=False)
    return b_splines


def b_splines_tmp(prec, u_expr, x_, y_, t_, t, degree, dt):
    '''
    :param degree: Degree of BSpline
    :param u_expr: analytic value of u
    :param dt: current time step size
    '''

    # grab keys
    keys = prec.read_data(0).keys()
    # use equidistant samples
    # apparently you need 2k+2 samples for degree k
    x_dist = float(dt) / (2 * degree + 2)
    nodes = []
    b_splines = {}
    # set times of samples
    for i in range(2*degree+3):
        if i == 0:
            nodes.append(0)
        elif i==2*degree+2:
            nodes.append(dt)
        else:
            nodes.append(x_dist * i)
    for k in keys:
        weights = []
        for i in range(2 * degree + 3):
            # get analytic solutions at the respective spatial and temporal location
            weights.append(u_expr.subs(x_, k[0]).subs(y_, k[1]).subs(t_,t+nodes[i]))
        # create b spline at (k[0], k[1])
        b1,b2,b3 = interpol.splrep(nodes, weights,s=0, k=degree)
        b_splines[k] = interpol.BSpline(b1, b2, b3, extrapolate=False)
    return b_splines
