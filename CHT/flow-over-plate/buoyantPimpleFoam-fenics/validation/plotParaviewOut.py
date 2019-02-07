from __future__ import division
import csv
from matplotlib import pyplot as plt
import numpy as np
from numpy import sqrt, exp
from scipy.special import erfc  # do not use erfc = 1 - numpy.erf, this leads to problems with floating point arithmetics!


T_inf = 300
T_c = 310
dT = T_c - T_inf

# symbols from [1] p.46, for geometry see [1] p.47
# parameters from [1] Fig.8 a)
Pr = 10 ** -2  # Prantl number
Re = 5 * 10 ** 2  # Reynolds number
b = 1  # length of slab
a = .25  # thickness of slab
lam = a / b  # aspect ratio
k = 1  # thermal conductivity ratio
tau = k / lam / sqrt(Re * Pr)  # dimensionless constant k/lam * sqrt(Re*Pr)

def read_and_plot(filename, label):

    reader = csv.DictReader(open(filename))

    result = {}
    for row in reader:
        print(row)
        for key in row:
            print(key)
            try:
                result[key].append(float(row[key]))
            except KeyError:
                result[key] = []
                result[key].append(float(row[key]))

    x = np.array(result['arc_length'])
    y = (np.array(result['T']) - T_inf) / dT
    plt.plot(x, y, label=label)


def plot_analytic(label):
    assert (Pr < 10)  # formula below only valid for Pr << 1; see [1] p.50
    assert (Re > 100)  # formula below only valid for Re >> 1; see [1] p.48
    theta_f = lambda X, Y: erfc(0.5 * Y / sqrt(X)) - exp(tau * (Y + tau * X)) * erfc(
        0.5 * Y / sqrt(X) + tau * sqrt(X))  # from [1] p. 51

    x = np.linspace(0, 1)
    T = theta_f(x,0)
    plt.plot(x, T, label=label)


read_and_plot('out_FE_OF.csv', "FE-OF coupling for Pr={Pr}, Re={Re}, k={k}".format(Pr=Pr, Re=Re, k=k))
read_and_plot('out_OF_OF.csv', "OF-OF coupling for Pr={Pr}, Re={Re}, k={k}".format(Pr=Pr, Re=Re, k=k))
plot_analytic("Analytical solution for Pr={Pr}, Re={Re}, k={k}".format(Pr=Pr, Re=Re, k=k))
plt.legend(loc='lower center')

plt.show()