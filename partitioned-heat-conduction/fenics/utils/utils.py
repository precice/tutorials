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
    nodes = np.linspace(0,dt,2*degree+3)  # <-- But these are 2k+3 samples
    weights = []
    b_splines = {}

    for node in nodes:
        weights.append(precice.read_data(node))

    for k in weights[0].keys():
        weights_k = []
        for i in range(2 * degree + 3):
            weights_k.append(weights[i][k])
        b1,b2,b3 = interpol.splrep(nodes, weights_k,s=0, k=degree)
        b_splines[k] = interpol.BSpline(b1, b2, b3, extrapolate=False)
    return b_splines
