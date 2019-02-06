# based on section 3 in [1]
# references:
# [1] Vynnycky, M., Kimura, S., Kanev, K., & Pop, I. (1998). Forced convection heat transfer from a flat plate: the conjugate problem. International Journal of Heat and Mass Transfer, 41(1), 45–59.

from __future__ import division
import numpy as np
from numpy import sqrt, exp
from scipy.special import erfc  # do not use erfc = 1 - numpy.erf, this leads to problems with floating point arithmetics!

# symbols from [1] p.46, for geometry see [1] p.47
# parameters fro [1] Fig.8 b)
Pr = 10**-2  # Prantl number
Re = 5*10**2  # Reynolds number
b = 1  # length of slab
a = .25  # thickness of slab
lam = a/b  # aspect ratio
k = 1  # thermal conductivity ratio
tau = k/lam/sqrt(Re*Pr)  # dimensionless constant k/lam * sqrt(Re*Pr)

assert (Pr < 10)  # formula below only valid for Pr << 1; see [1] p.50
assert (Re > 100)  # formula below only valid for Re >> 1; see [1] p.48
theta_f = lambda X, Y: erfc(0.5 * Y / sqrt(X)) - exp(tau * (Y + tau * X)) * erfc(0.5 * Y / sqrt(X) + tau * sqrt(X))  # from [1] p. 51
theta_bb = lambda X: 1 - exp(tau**2 * X) * erfc(tau * sqrt(X))  # from [1] p. 51

x = np.linspace(-.5, .5, 100)
X = x + .5

from collections import defaultdict

# we use a nested dict as data structure (https://stackoverflow.com/a/8702435/5158031)
nested_dict = lambda: defaultdict(nested_dict)
vynnycky_analytical = nested_dict()  # we put data in this dict. Vynnycky[Pr][Re][k][coordinate]
# Fig 8 a)
_Pr, _Re, _k = 10**-2, 5*10**2, 1
vynnycky_analytical[_Pr][_Re][_k]['x'] =     np.array([0  , 27 , 53 , 136, 202, 276, 316, 343, 406, 457, 468, 484, 494]) / 498 - .5
vynnycky_analytical[_Pr][_Re][_k]['y'] = 1 - np.array([170, 160, 148, 125, 113, 104, 100, 99 , 97 , 98 , 100, 102, 103]) / 310
_Pr, _Re, _k = 10**-2, 5*10**2, 5
vynnycky_analytical[_Pr][_Re][_k]['x'] =     np.array([11, 42, 73, 102, 132, 163, 193, 223, 254, 285, 314, 345, 376, 405, 436, 467, 497]) /498 - .5
vynnycky_analytical[_Pr][_Re][_k]['y'] = 1 - np.array([79, 65, 55, 49 , 44 , 40 , 35 , 33 , 31 , 29 , 27 , 25 , 24 , 23 , 22 , 22 , 22]) / 310
_Pr, _Re, _k = 10**-2, 5*10**2, 20
vynnycky_analytical[_Pr][_Re][_k]['x'] =     np.array([11, 42, 73, 102, 132, 163, 193, 223, 254, 285, 314, 345, 376, 405, 436, 467, 497]) /498 - .5
vynnycky_analytical[_Pr][_Re][_k]['y'] = 1 - np.array([24, 20, 17, 15 , 13 , 12 , 10 , 9  , 8  , 7  , 7  , 6.5, 6  , 5.5, 5  , 5  , 5  ]) / 310
# Fig 8 b)
_Pr, _Re, _k = 10**-2, 10**4, 20
vynnycky_analytical[_Pr][_Re][_k]['x'] = np.array([21, 89, 155, 220, 285, 349, 416, 480, 546, 613, 679, 743, 810, 875, 940, 1007, 1071]) / 1071 - .5
vynnycky_analytical[_Pr][_Re][_k]['y'] = 1 - np.array([204, 164, 140, 124, 110, 98, 90, 82, 77, 71, 68, 64, 61, 58, 58, 58, 55]) / 752

theta_b = theta_bb(X)

from matplotlib import pyplot as plt
from scipy.signal import savgol_filter

plt.plot(x,theta_b, label=r"analytical solution $\theta_b(x)$")
y = vynnycky_analytical[Pr][Re][k]['y']
plt.plot(vynnycky_analytical[Pr][Re][k]['x'], savgol_filter(y, 11, 3), 'o',
         label="Vynnycky, analytical for Pr={Pr}, Re={Re}, k={k}".format(Pr=Pr, Re=Re, k=k))
plt.plot(vynnycky_analytical[Pr][Re][k]['x'], y, 'sb',
         label="Vynnycky, analytical for Pr={Pr}, Re={Re}, k={k}".format(Pr=Pr, Re=Re, k=k))
plt.xlabel("x")
plt.ylabel(r"$\theta_b$")
plt.legend(loc='lower center')
plt.show()

""" Fig.2 & Fig.3
x, y = np.meshgrid(np.linspace(-.5,.5,100), np.linspace(0,3,100))
X = x + .5
Y = y + 0

val = theta_f(X,Y)

plt.contour(x,y,val)
plt.show()
"""
