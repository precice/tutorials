from typing import Tuple, List
import numpy as np
import numbers
import scipy as sp
from enum import Enum


class TimeSteppingSchemes(Enum):
    NEWMARK_BETA = "Newmark_beta"
    GENERALIZED_ALPHA = "generalized_alpha"
    RUNGE_KUTTA_4 = "runge_kutta_4"
    Radau_IIA = "radauIIA"  # 5th order implicit


class GeneralizedAlpha():
    alpha_f = None
    alpha_m = None
    gamma = None
    beta = None
    mass = None
    stiffness = None

    def __init__(self, stiffness, mass, alpha_f=0.4, alpha_m=0.2) -> None:
        self.alpha_f = alpha_f
        self.alpha_m = alpha_m

        self.gamma = 0.5 - self.alpha_m + self.alpha_f
        self.beta = 0.25 * (self.gamma + 0.5)

        self.stiffness = stiffness
        self.mass = mass

    def do_step(self, u0, v0, a0, rhs, dt) -> Tuple[float, float, float]:
        f = rhs((1 - self.alpha_f) * dt)

        m = 3 * [None]
        m[0] = (1 - self.alpha_m) / (self.beta * dt**2)
        m[1] = (1 - self.alpha_m) / (self.beta * dt)
        m[2] = (1 - self.alpha_m - 2 * self.beta) / (2 * self.beta)

        k_bar = self.stiffness * (1 - self.alpha_f) + m[0] * self.mass

        # do generalized alpha step
        if (type(self.stiffness)) is np.ndarray:
            u1 = np.linalg.solve(
                k_bar,
                (f - self.alpha_f * self.stiffness.dot(u0) + self.mass.dot((m[0] * u0 + m[1] * v0 + m[2] * a0)))
            )
        else:
            u1 = (f - self.alpha_f * self.stiffness * u0 + self.mass * (m[0] * u0 + m[1] * v0 + m[2] * a0)) / k_bar

        a1 = 1.0 / (self.beta * dt**2) * (u1 - u0 - dt * v0) - (1 - 2 * self.beta) / (2 * self.beta) * a0
        v1 = v0 + dt * ((1 - self.gamma) * a0 + self.gamma * a1)

        return u1, v1, a1


class RungeKutta4():
    a = np.array([[0, 0, 0, 0],
                  [0.5, 0, 0, 0],
                  [0, 0.5, 0, 0],
                  [0, 0, 1.0, 0]])
    b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])
    c = np.array([0, 0.5, 0.5, 1])

    def __init__(self, ode_system) -> None:
        self.ode_system = ode_system
        pass

    def do_step(self, u0, v0, _, rhs, dt) -> Tuple[float, float, float]:
        assert (isinstance(u0, type(v0)))

        k = 4 * [None]  # store stages in list

        if isinstance(u0, np.ndarray):
            x0 = np.concatenate([u0, v0])
            def f(t,x): return self.ode_system.dot(x) + np.concatenate([np.array([0, 0]), rhs(t)])
        elif isinstance(u0, numbers.Number):
            x0 = np.array([u0, v0])
            def f(t,x): return self.ode_system.dot(x) + np.array([0, rhs(t)])
        else:
            raise Exception(f"Cannot handle input type {type(u0)} of u and v")

        k[0] = f(self.c[0] * dt, x0)
        k[1] = f(self.c[1] * dt, x0 + self.a[1, 0] * k[0] * dt)
        k[2] = f(self.c[2] * dt, x0 + self.a[2, 1] * k[1] * dt)
        k[3] = f(self.c[3] * dt, x0 + self.a[3, 2] * k[2] * dt)

        x1 = x0 + dt * sum(b_i * k_i for k_i, b_i in zip(k, self.b))

        if isinstance(u0, np.ndarray):
            u1 = x1[0:2]
            v1 = x1[2:4]
        elif isinstance(u0, numbers.Number):
            u1 = x1[0]
            v1 = x1[1]

        return u1, v1, None


class RadauIIA():
    def __init__(self, ode_system) -> None:
        self.ode_system = ode_system
        pass

    def do_step(self, u0, v0, _, rhs, dt) -> Tuple[float, float, float]:
        assert (isinstance(u0, type(v0)))

        t0 = 0

        if isinstance(u0, np.ndarray):
            x0 = np.concatenate([u0, v0])
            def f(t,x): return self.ode_system.dot(x) + np.concatenate([np.array([np.zeros_like(t), np.zeros_like(t)]), rhs(t)])
        elif isinstance(u0, numbers.Number):
            x0 = np.array([u0, v0])
            def f(t,x): return self.ode_system.dot(x) + np.array([np.zeros_like(t), rhs(t)])
        else:
            raise Exception(f"Cannot handle input type {type(u0)} of u and v")

        # use adaptive time stepping; dense_output=True allows us to sample from continuous function later
        ret = sp.integrate.solve_ivp(f, [t0, t0 + dt], x0, method="Radau",
                                     dense_output=True, rtol=10e-5, atol=10e-9)

        if isinstance(u0, np.ndarray):
            u1, v1 = ret.y[0:2, -1], ret.y[2:4, -1]
        elif isinstance(u0, numbers.Number):
            u1, v1 = ret.y[:, -1]

        return u1, v1, None, ret.sol
