# based on https://github.com/precice/elastictube1d and
# [1] J. Degroote, P. Bruggeman, R. Haelterman, and J. Vierendeels. Stability of a coupling technique for partitioned solvers in FSI applications. Computers & Structures, 2008.
# for time integration details see
# [2] Gresho, P. M., & Sani, R. L. (2000). Incompressible Flow and the Finite Element Method, Isothermal Laminar Flow. John Wiley & Sons. Retrieved from http://books.google.de/books?id=m_tQAAAAMAAJ

from __future__ import division, print_function
import numpy as np


def perform_partitioned_theta_scheme_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in,
                                          custom_coupling, theta=1):

    k = 0

    # initial guess for Newtons method
    pressure1 = np.copy(pressure0)
    velocity1 = np.copy(velocity0)

    crossSection_couple = 2 * [None]
    if custom_coupling:
        # set cross sections corresponding to point in time
        crossSection_couple[0] = crossSection0
        crossSection_couple[1] = crossSection1
    else:
        # set both cross sections equal to input -> depending on input: implicit or explicit coupling
        crossSection_couple[0] = crossSection1
        crossSection_couple[1] = crossSection1

    N = pressure0.shape[0] - 1

    alpha = 0
    success = True

    E = 10000  # elasticity module
    c_mk = np.sqrt(E / 2 * np.sqrt(np.pi))  # wave speed

    while success:  # perform Newton iterations to solve nonlinear system of equations

        # compute residual
        res = np.zeros(2 * N + 2)

        for i in range(1, N):
            # Momentum
            res[i] = (velocity0[i] * crossSection0[i] -
                      velocity1[i] * crossSection1[i]) * dx / tau

            res[i] += .25 * theta * (- crossSection_couple[1][i + 1] * velocity1[i] * velocity1[i + 1]
                                     - crossSection_couple[1][i] * velocity1[i] * velocity1[i + 1])
            res[i] += .25 * (1 - theta) * (- crossSection_couple[0][i + 1] * velocity0[i] * velocity0[i + 1]
                                           - crossSection_couple[0][i] * velocity0[i] * velocity0[i + 1])

            res[i] += .25 * theta * (- crossSection_couple[1][i + 1] * velocity1[i] * velocity1[i]
                                     - crossSection_couple[1][i] * velocity1[i] * velocity1[i]
                                     + crossSection_couple[1][i] * velocity1[i - 1] * velocity1[i]
                                     + crossSection_couple[1][i - 1] * velocity1[i - 1] * velocity1[i])
            res[i] += .25 * (1 - theta) * (- crossSection_couple[0][i + 1] * velocity0[i] * velocity0[i]
                                           - crossSection_couple[0][i] * velocity0[i] * velocity0[i]
                                           + crossSection_couple[0][i] * velocity0[i - 1] * velocity0[i]
                                           + crossSection_couple[0][i - 1] * velocity0[i - 1] * velocity0[i])

            res[i] += .25 * theta * (+ crossSection_couple[1][i - 1] * velocity1[i - 1] * velocity1[i - 1]
                                     + crossSection_couple[1][i] * velocity1[i - 1] * velocity1[i - 1])
            res[i] += .25 * (1 - theta) * (+ crossSection_couple[0][i - 1] * velocity0[i - 1] * velocity0[i - 1]
                                           + crossSection_couple[0][i] * velocity0[i - 1] * velocity0[i - 1])

            res[i] += .25 * theta * (+ crossSection_couple[1][i - 1] * pressure1[i - 1]
                                     + crossSection_couple[1][i] * pressure1[i - 1]
                                     - crossSection_couple[1][i - 1] * pressure1[i]
                                     + crossSection_couple[1][i + 1] * pressure1[i]
                                     - crossSection_couple[1][i] * pressure1[i + 1]
                                     - crossSection_couple[1][i + 1] * pressure1[i + 1])
            res[i] += .25 * (1 - theta) * (+ crossSection_couple[0][i - 1] * pressure0[i - 1]
                                           + crossSection_couple[0][i] * pressure0[i - 1]
                                           - crossSection_couple[0][i - 1] * pressure0[i]
                                           + crossSection_couple[0][i + 1] * pressure0[i]
                                           - crossSection_couple[0][i] * pressure0[i + 1]
                                           - crossSection_couple[0][i + 1] * pressure0[i + 1])

            # Continuity (we only care about values at n+1, see [2],p.737,eq.(3.16-25))
            res[i + N + 1] = (crossSection0[i] - crossSection1[i]) * dx / tau
            res[i + N + 1] += .25 * theta * (+ crossSection_couple[1][i - 1] * velocity1[i - 1]
                                             + crossSection_couple[1][i] * velocity1[i - 1]
                                             + crossSection_couple[1][i - 1] * velocity1[i]
                                             - crossSection_couple[1][i + 1] * velocity1[i]
                                             - crossSection_couple[1][i] * velocity1[i + 1]
                                             - crossSection_couple[1][i + 1] * velocity1[i + 1])
            res[i + N + 1] += .25 * (1 - theta) * (+ crossSection_couple[0][i - 1] * velocity0[i - 1]
                                                   + crossSection_couple[0][i] * velocity0[i - 1]
                                                   + crossSection_couple[0][i - 1] * velocity0[i]
                                                   - crossSection_couple[0][i + 1] * velocity0[i]
                                                   - crossSection_couple[0][i] * velocity0[i + 1]
                                                   - crossSection_couple[0][i + 1] * velocity0[i + 1])
            res[i + N + 1] += alpha * theta * (pressure1[i - 1] - 2 * pressure1[i] + pressure1[i + 1])

        # Boundary

        # Velocity Inlet is prescribed
        res[0] = velocity_in - velocity1[0]

        # Pressure Inlet is lineary interpolated
        res[N + 1] = -pressure1[0] + 2 * pressure1[1] - pressure1[2]

        # Velocity Outlet is lineary interpolated
        res[N] = -velocity1[-1] + 2 * velocity1[-2] - velocity1[-3]

        # Pressure Outlet is "non-reflecting"
        tmp2 = np.sqrt(c_mk ** 2 - pressure0[-1] / 2) - (velocity1[-1] - velocity0[-1]) / 4
        res[2 * N + 1] = -pressure1[-1] + 2 * (c_mk ** 2 - tmp2 * tmp2)

        k += 1  # Iteration Count

        # compute relative norm of residual
        norm_1 = np.sqrt(res.dot(res))
        norm_2 = np.sqrt(pressure1.dot(pressure1) + velocity1.dot(velocity1))
        norm = norm_1 / norm_2

        if norm < 1e-10 and k > 1:
            break  # Nonlinear Solver success
        elif k > 1000:
            print(
                "Nonlinear Solver break, iterations: %i, residual norm: %e\n" % (k, norm))
            velocity1[:] = np.nan
            pressure1[:] = np.nan
            success = False
            break
        # else:
        # perform another iteration of newton's method

        # compute Jacobian for Newton's method
        system = np.zeros([N + N + 2, N + N + 2])

        for i in range(1, N):
            # Momentum, Velocity see [1] eq. (13b)
            system[i][i - 1] += .25 * theta * (- 2 * crossSection_couple[1][i - 1] * velocity1[i - 1]
                                               - 2 * crossSection_couple[1][i] * velocity1[i - 1]
                                               - crossSection_couple[1][i] * velocity1[i]
                                               - crossSection_couple[1][i - 1] * velocity1[i])
            system[i][i] += crossSection1[i] * dx / tau
            system[i][i] += .25 * theta * (+ crossSection_couple[1][i + 1] * velocity1[i + 1]
                                           + crossSection_couple[1][i] * velocity1[i + 1]
                                           + crossSection_couple[1][i + 1] * velocity1[i] * 2
                                           + crossSection_couple[1][i] * velocity1[i] * 2
                                           - crossSection_couple[1][i] * velocity1[i - 1]
                                           - crossSection_couple[1][i - 1] * velocity1[i - 1])
            system[i][i + 1] += .25 * theta * (crossSection_couple[1][i + 1] * velocity1[i]
                                               + crossSection_couple[1][i] * velocity1[i])

            # Momentum, Pressure see [1] eq. (13b)
            system[i][N + 1 + i - 1] += .25 * theta * (- crossSection_couple[1][i - 1] - crossSection_couple[1][i])
            system[i][N + 1 + i] += .25 * theta * (+ crossSection_couple[1][i - 1] - crossSection_couple[1][i + 1])
            system[i][N + 1 + i + 1] += .25 * theta * (+ crossSection_couple[1][i] + crossSection_couple[1][i + 1])

            # Continuity, Velocity see [1] eq. (13a)
            system[i + N + 1][i - 1] += .25 * theta * (- crossSection_couple[1][i - 1] - crossSection_couple[1][i])
            system[i + N + 1][i] += .25 * theta * (- crossSection_couple[1][i - 1] + crossSection_couple[1][i + 1])
            system[i + N + 1][i + 1] += .25 * theta * (+ crossSection_couple[1][i] + crossSection_couple[1][i + 1])

            # Continuity, Pressure see [1] eq. (13a)
            system[i + N + 1][N + 1 + i - 1] += - alpha * theta
            system[i + N + 1][N + 1 + i] += 2 * alpha * theta
            system[i + N + 1][N + 1 + i + 1] += - alpha * theta

        # Velocity Inlet is prescribed
        system[0][0] = 1
        # Pressure Inlet is linearly interpolated [1] eq. (14a)
        system[N + 1][N + 1] = 1
        system[N + 1][N + 2] = -2
        system[N + 1][N + 3] = 1
        # Velocity Outlet is linearly interpolated [1] eq. (14b)
        system[N][N] = 1
        system[N][N - 1] = -2
        system[N][N - 2] = 1

        # Pressure Outlet is Non-Reflecting [1] eq. (15)
        system[2 * N + 1][2 * N + 1] = 1
        system[2 * N + 1][N] = -(np.sqrt(c_mk ** 2 - pressure0[-1] / 2) - (velocity1[-1] - velocity0[-1]) / 4)

        try:
            solution = np.linalg.solve(system, res)
        except np.linalg.LinAlgError:
            print("LINALGERROR! SINGULAR MATRIX")
            velocity1[:] = np.nan
            pressure1[:] = np.nan
            success = False
            break

        velocity1 += solution[:N + 1]
        pressure1 += solution[N + 1:]

    return velocity1, pressure1, success


def perform_partitioned_implicit_euler_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau,
                                            velocity_in, custom_coupling):
    return perform_partitioned_theta_scheme_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau,
                                                 velocity_in, custom_coupling=False, theta=1)


def perform_partitioned_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau,
                                                       velocity_in, custom_coupling):
    return perform_partitioned_theta_scheme_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau,
                                                 velocity_in, custom_coupling, theta=.5)
