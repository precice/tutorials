from typing import Tuple, Callable
import numpy as np
import scipy as sp
from scipy.integrate import OdeSolution
from enum import Enum


class TimeSteppingSchemes(Enum):
    NEWMARK_BETA = "Newmark_beta"
    GENERALIZED_ALPHA = "generalized_alpha"
    RUNGE_KUTTA_4 = "runge_kutta_4"
    Radau_IIA = "radauIIA"  # 5th order implicit


class TimeStepper:
    """Time stepper for a mass spring system with given mass and spring stiffness
    """

    def __init__(self,
                 stiffness: float,
                 mass: float) -> None:
        """Initializes time stepper with parameters describing mass spring system

        Args:
            stiffness (float): the stiffness of the spring connected to the mass
            mass (float): the mass
        """
        raise NotImplementedError()

    def do_step(self,
                u0: float,
                v0: float,
                a0: float,
                rhs: Callable[[float], float],
                dt: float
                ) -> Tuple[float, float, float]:
        """Perform a time step of size dt with given time stepper

        Args:
            u0 (float): displacement at time t0
            v0 (float): velocity at time t0
            a0 (float): acceleration at time t0
            rhs (Callable[[float],float]): time dependent right-hand side
            dt (float): time step size

        Returns:
            Tuple[float, float, float]: returns computed displacement, velocity, and acceleration at time t1 = t0 + dt.
        """
        raise NotImplementedError()


class GeneralizedAlpha(TimeStepper):
    """TimeStepper implementing generalized Alpha or Newmark Beta scheme (depends on parameters alpha_f and alpha_m set in constructor)
    """
    alpha_f: float
    alpha_m: float

    def __init__(self, stiffness: float, mass: float, alpha_f: float = 0.4, alpha_m: float = 0.2) -> None:
        self.alpha_f = alpha_f
        self.alpha_m = alpha_m

        self.gamma = 0.5 - self.alpha_m + self.alpha_f
        self.beta = 0.25 * (self.gamma + 0.5)

        self.stiffness = stiffness
        self.mass = mass

    def do_step(self,
                u0: float,
                v0: float,
                a0: float,
                rhs: Callable[[float], float],
                dt: float
                ) -> Tuple[float, float, float]:
        f = rhs((1.0 - self.alpha_f) * dt)

        m = [(1 - self.alpha_m) / (self.beta * dt**2),
             (1 - self.alpha_m) / (self.beta * dt),
             (1 - self.alpha_m - 2 * self.beta) / (2 * self.beta)]

        k_bar = self.stiffness * (1 - self.alpha_f) + m[0] * self.mass

        # do generalized alpha step
        u1 = (f - self.alpha_f * self.stiffness * u0 + self.mass * (m[0] * u0 + m[1] * v0 + m[2] * a0)) / k_bar
        a1 = 1.0 / (self.beta * dt**2) * (u1 - u0 - dt * v0) - (1 - 2 * self.beta) / (2 * self.beta) * a0
        v1 = v0 + dt * ((1 - self.gamma) * a0 + self.gamma * a1)

        return u1, v1, a1


class RungeKutta4(TimeStepper):
    """TimeStepper implementing classic Runge Kutta scheme (4 stages)
    """
    # parameters from Butcher tableau of classic Runge Kutta scheme (RK4)
    a = np.array([[0, 0, 0, 0],
                  [0.5, 0, 0, 0],
                  [0, 0.5, 0, 0],
                  [0, 0, 1.0, 0]])
    b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])
    c = np.array([0, 0.5, 0.5, 1])

    def __init__(self, stiffness: float, mass: float) -> None:
        self.ode_system = np.array([
            [0, 1],  # du
            [-stiffness / mass, 0],  # dv
        ])

    def do_step(self,
                u0: float,
                v0: float,
                _: float,
                rhs: Callable[[float], float],
                dt: float
                ) -> Tuple[float, float, float]:
        k = 4 * [None]  # store stages in list

        x0 = np.array([u0, v0])
        def f(t, x): return self.ode_system.dot(x) + np.array([0, rhs(t)])

        k[0] = f(self.c[0] * dt, x0)
        k[1] = f(self.c[1] * dt, x0 + self.a[1, 0] * k[0] * dt)
        k[2] = f(self.c[2] * dt, x0 + self.a[2, 1] * k[1] * dt)
        k[3] = f(self.c[3] * dt, x0 + self.a[3, 2] * k[2] * dt)

        x1 = x0 + dt * sum(b_i * k_i for k_i, b_i in zip(k, self.b))

        return x1[0], x1[1], f(dt, x1)[1]


class RadauIIA(TimeStepper):
    """Perform a step with the adaptive RadauIIA time stepper of scipy.integrate
    """
    _dense_output: OdeSolution

    def __init__(self, stiffness: float, mass: float) -> None:
        self.ode_system = np.array([
            [0, 1],  # du
            [-stiffness / mass, 0],  # dv
        ])

    def do_step(self,
                u0: float,
                v0: float,
                _: float,
                rhs: Callable[[float], float],
                dt: float
                ) -> Tuple[float, float, float]:
        t0, t1 = 0, dt

        x0 = np.array([u0, v0])
        def f(t, x): return self.ode_system.dot(x) + np.array([np.zeros_like(t), rhs(t)])

        # use adaptive time stepping; dense_output=True allows us to sample from continuous function later
        ret = sp.integrate.solve_ivp(f, [t0, t1], x0, method="Radau",
                                     dense_output=True, rtol=10e-5, atol=10e-9)

        self._dense_output = ret.sol  # store dense output in class

        return ret.y[0, -1], ret.y[1, -1], f(dt, ret.y[:, -1])[1]

    @property
    def dense_output(self) -> OdeSolution:
        """Returns dense output created for the last call of do_step

        Returns:
            OdeSolution: dense output
        """
        return self._dense_output
