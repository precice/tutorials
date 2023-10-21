import argparse
import sys
import numpy as np

from fenics import dx, TrialFunction, TestFunction, solve, inner, grad

def approx_derivative(precice, x, max_x):
    """
    :param precice: precice object of the problem
    :param x: where the derivative should be approximated. 0<=x<=max_x
    :param max_x: the maximum x which is supported by the precice interpolant (it's precice-dt!)
    approximates the derivative at time x
    """
    if x < 0 or max_x < x:
        raise Exception("x is out of bounds")
        return {}
    elif x == 0:
        return approx_derivative_forw(precice, x, max_x)
    elif x == max_x:
        return approx_derivative_backw(precice, x, max_x)

    derivatives = {}
    # convergence order of 16
    order = 8
    # minimal difference between d_n and d_np1 to iterate again
    delta_tol = 1e-15
    # last difference between d_n and d_np1
    delta_last = sys.float_info.max
    # last derivative approximation
    d_n = sys.float_info.max
    d_np1 = 0
    # initial step size
    h = 0.5
    # iterate over all cells on the boundary
    for k in precice.read_data(x).keys():
        # reset variables for stopping criterion
        delta_last = sys.float_info.max
        d_n = sys.float_info.max
        d_np1 = 0

        # initialize system of equations
        # Todo: for each step size h, the matrix is for all evaluations at k
        #  identical. It can be build once and then used with QR decomposition
        dim = order * 2 + 1
        b = [0] * dim
        A = [[0] * dim for i in range(dim)]
        factorials = []
        for i in range(dim):
            factorials.append(np.math.factorial(i))

        while delta_tol <= abs(d_np1 - d_n):
            # build matrix
            for i in range(dim):
                for j in range(dim):
                    A[i][j] = pow((i - order) * h, j) / factorials[j]
            # build b
            # flag is set to false, if at least one precice call is out of bounds, i.e. not in [0, max_x]
            # the approximation will then continue with h/2 and not be stopped
            flag = True
            for i in range(dim):
                if 0 <= x + (i - order) * h <= max_x:
                    # get the value to the position k at time x + (i - order) * h
                    b[i] = precice.read_data(x + (i - order) * h)[k]
                else:
                    flag = False
                    break
            if flag:
                # solve lse
                tmp = np.linalg.solve(A, b)
                # check if delta is greater than the last iteration
                if delta_last < abs(d_np1 - tmp[1]):
                    # delta is greater than the last iteration
                    # -> when doing further iterations the approx. will get worse
                    derivatives[k] = d_np1
                    break
                else:
                    delta_last = abs(d_np1 - tmp[1])

                # update values
                d_n = d_np1
                d_np1 = tmp[1]
            h /= 2
            # avoid that h gets too small to evade instabilities
            if h < 1e-12:
                # just use the current approximation if it gets to small
                derivatives[k] = d_np1
                break
        derivatives[k] = d_np1
    return derivatives


def approx_derivative_forw(precice, x, max_x):
    """
    Same functionality as approx_derivative. Difference is, that forward difference quotient is used
    instead of central difference quotient. If avoidable use the other function!
    :param precice: precice object of the problem
    :param x: where the derivative should be approximated. 0<=x<=max_x
    :param max_x: the maximum x which is supported by the precice interpolant (it's precice-dt!)
    """
    # how many derivatives are needed
    length = len(precice.read_data(x))
    if x < 0 or max_x <= x:
        raise Exception("x is out of bounds")
        return {}

    derivatives = {}
    # convergence order of 8
    order = 8
    # minimal difference between d_n and d_np1 to iterate again
    delta_tol = 1e-15
    # last difference between d_n and d_np1
    delta_last = sys.float_info.max
    # last derivative approximation
    d_n = sys.float_info.max
    d_np1 = 0
    # initial step size
    h = 0.5
    # iterate over all cells on the boundary
    for k in precice.read_data(x).keys():
        # reset variables for stopping criterion
        delta_last = sys.float_info.max
        d_n = sys.float_info.max
        d_np1 = 0

        # initialize system of equations
        # Todo: for each step size h, the matrix is for all evaluations at k
        #  identical. It can be build once and then used with QR decomposition
        dim = order * 2 + 1
        b = [0] * dim
        A = [[0] * dim for i in range(dim)]
        factorials = []
        for i in range(dim):
            factorials.append(np.math.factorial(i))

        while delta_tol <= abs(d_np1 - d_n):
            # build matrix
            for i in range(dim):
                for j in range(dim):
                    A[i][j] = pow(i * h, j) / factorials[j]
            # build b
            # flag is set to false, if at least one precice call is out of bounds, i.e. not in [0, max_x]
            # the approximation will then continue with h/2 and not be stopped
            flag = True
            for i in range(dim):
                if x + i * h <= max_x:
                    # get the value to the position k at time x + i * h
                    # as we use forward difference quotient, the sampled values are only >=x
                    b[i] = precice.read_data(x + i * h)[k]
                else:
                    flag = False
                    break
            if flag:
                # solve lse
                tmp = np.linalg.solve(A, b)
                # check if delta is greater than the last iteration
                if delta_last < abs(d_np1 - tmp[1]):
                    # delta is greater than the last iteration
                    # -> when doing further iterations the approx. will get worse
                    derivatives[k] = d_np1
                    break
                else:
                    delta_last = abs(d_np1 - tmp[1])

                # update values
                d_n = d_np1
                d_np1 = tmp[1]
            h /= 2
            # avoid that h gets too small to evade instabilities
            if h < 1e-12:
                # just use the current approximation if it gets to small
                derivatives[k] = d_np1
                break
        derivatives[k] = d_np1
    return derivatives


def approx_derivative_backw(precice, x, max_x):
    """
    Same functionality as approx_derivative. Difference is, that backward difference quotient is used
    instead of central difference quotient. If avoidable use the other function!
    :param precice: precice object of the problem
    :param x: where the derivative should be approximated. 0<=x<=max_x
    :param max_x: the maximum x which is supported by the precice interpolant (it's precice-dt!)
    """
    # how many derivatives are needed
    length = len(precice.read_data(x))
    if x <= 0 or max_x < x:
        raise Exception("x is out of bounds")
        return {}

    derivatives = {}
    # convergence order of 8
    order = 8
    # minimal difference between d_n and d_np1 to iterate again
    delta_tol = 1e-15
    # last difference between d_n and d_np1
    delta_last = sys.float_info.max
    # last derivative approximation
    d_n = sys.float_info.max
    d_np1 = 0
    # initial step size
    h = 0.5
    # iterate over all cells on the boundary
    for k in precice.read_data(x).keys():
        # reset variables for stopping criterion
        h = 0.5
        delta_last = sys.float_info.max
        d_n = sys.float_info.max
        d_np1 = 0

        # initialize system of equations
        # Todo: for each step size h, the matrix is for all evaluations at k
        #  identical. It can be build once and then used with QR decomposition
        dim = order * 2 + 1
        b = [0] * dim
        A = [[0] * dim for i in range(dim)]
        factorials = []
        for i in range(dim):
            factorials.append(np.math.factorial(i))

        while delta_tol <= abs(d_np1 - d_n):
            # build matrix
            for i in range(dim):
                for j in range(dim):
                    A[i][j] = pow(-i * h, j) / factorials[j]
            # build b
            # flag is set to false, if at least one precice call is out of bounds, i.e. not in [0, max_x]
            # the approximation will then continue with h/2 and not be stopped
            flag = True
            for i in range(dim):
                if x - i * h >= 0:
                    # get the value to the position k at time x - i * h
                    # as we use forward difference quotient, the sampled values are only <=x
                    b[i] = precice.read_data(x - i * h)[k]
                else:
                    flag = False
                    break
            if flag:
                # solve lse
                tmp = np.linalg.solve(A, b)
                # check if delta is greater than the last iteration
                if delta_last < abs(d_np1 - tmp[1]):
                    # delta is greater than the last iteration
                    # -> when doing further iterations the approx. will get worse
                    derivatives[k] = d_np1
                    break
                else:
                    delta_last = abs(d_np1 - tmp[1])

                # update values
                d_n = d_np1
                d_np1 = tmp[1]
            h /= 2
            # avoid that h gets too small to evade instabilities
            if h < 1e-12:
                # just use the current approximation if it gets to small
                derivatives[k] = d_np1
                break
        derivatives[k] = d_np1
    return derivatives
