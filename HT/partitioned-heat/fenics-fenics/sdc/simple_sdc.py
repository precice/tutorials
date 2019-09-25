from sympy.integrals.quadrature import gauss_lobatto
import numpy as np


# Generate nodes and weights
M = 3
x, w = gauss_lobatto(M, n_digits=30)

# number of correction sweeps
K = 4

# expected order
p = 2*M-2

# convert x to numpy array
x = np.array([float(xx) for xx in x])
w = np.array([float(ww) for ww in w])

# Transform to [0, 1]
x = 0.5*(x + 1.0)
w *= 0.5

Q = np.zeros((M, M))
for j in range(M):
    y = np.zeros(M)
    y[j] = 1.0
    # c is now the j-th Lagrange polynomial
    c = np.polyfit(np.float64(x), y, M-1)
    # since c is a polynomial, we can integrate it exactly
    cint = np.polyint(c)
    for m in range(M):
        Q[m, j] = np.polyval(cint, x[m]) - np.polyval(cint, 0.0)


def sdc_step(y0, t0, black_box_implicit_euler, f, dt, V):
    from fenics import Function
    # initialize with implicit euler
    y = [[Function(V) for _ in range(M)] for _ in range(K+1)]
    for i in range(K+1):
        y[i][0].assign(y0)
    d_tau = (x[1:] - x[:-1]) * dt
    t = t0 + x * dt
    k = 0
    for mm in range(M-1):
        y[k][mm+1].assign(black_box_implicit_euler(y[k][mm], t[mm], d_tau[mm], mm))

    S = Q[1:][:] - Q[:-1][:]
    S *= dt
    for k in range(K):
        buffered_f = []
        for i in range(M):
            buffered_f.append(f(y[k][i], t[i], i))
        for mm in range(M-1):
            correction_sum = 0
            for i in range(M):
                u_rhs = Function(V)
                u_rhs.assign(S[mm][i] * buffered_f[i])  # = S[mm][i] * f(y[k][i], t[i])
                correction_sum += u_rhs
            u_rhs_last = f(y[k][mm+1], t[mm+1], mm+1)
            corrected_y = Function(V)
            corrected_y.assign(y[k+1][mm] - d_tau[mm] * u_rhs_last + correction_sum)
            y[k + 1][mm + 1].assign(black_box_implicit_euler(corrected_y,
                                                             t[mm],
                                                             d_tau[mm],
                                                             mm))

    return y[-1][-1]
