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
        nodes.append(x_dist * i)
        weights.append(precice.read_data(x_dist * i))
    for k in weights[0].keys():
        for i in range(2 * degree + 3):
            weights_k = []
            weights_k.append(weights[i][k])
        b_splines[k] = interpol.BSpline(nodes, weights_k, degree)
    return b_splines
