import os
from coupling_schemes import CouplingScheme
import itertools
from tabulate import tabulate
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np

def create_error_table(prefix, evaluated_wr, evaluated_dT, coupling_schemes):

    simulationTime = 10.0
    data_D = []
    data_D.append(["WR", "dT", "cpl", "error"])
    data_N = []
    data_N.append(["WR", "dT", "cpl", "error"])

    for wr, dT, coupling_scheme in itertools.product(evaluated_wr, evaluated_dT, coupling_schemes):
        wr_tag = "WR{wr}".format(wr=wr)
        window_tag = "dT{dT}".format(dT=dT)
        folder = os.path.join(prefix, wr_tag, window_tag, coupling_scheme)
        path_D = os.path.join(folder,"errors*DIRICHLET")
        path_N = os.path.join(folder,"errors*NEUMANN")
        assert(glob.glob(path_D).__len__() <= 1)
        if(glob.glob(path_D).__len__() == 1):
            try:
                with open(glob.glob(path_D)[0]) as file:
                    for line in file.readlines():
                        lastline = line
                    lastline = lastline.split(", ")
                    dt = float(lastline[0])
                    error = float(lastline[1])
                    data_D.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "error": error})
            except FileNotFoundError:
                data_D.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "error": "-"})
        else:
            data_D.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "error": "-"})

        assert(glob.glob(path_N).__len__() <= 1)
        if(glob.glob(path_N).__len__() == 1):
            try:
                with open(glob.glob(path_N)[0]) as file:
                    for line in file.readlines():
                        lastline = line
                    lastline = lastline.split(", ")
                    dt = float(lastline[0])
                    error = float(lastline[1])
                    data_N.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "error": error})
            except FileNotFoundError:
                data_N.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "error": "-"})
        else:
            data_N.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "error": "-"})

    table = []
    keys = ["dT{dT}".format(dT=dT) for dT in evaluated_dT]
    table.append(["error"] + keys)

    structured_data = dict()

    for d_D, d_N in zip(data_D[1:],data_N[1:]):
        wrcplkey = "WR{wr}, cpl={cpl}".format(wr=d_D["WR"], cpl=d_D["cpl"])

        try:
            structured_data[wrcplkey]
        except KeyError:
            structured_data[wrcplkey] = dict()
        dTkey = "dT{dT}".format(dT=d_D["dT"])
        try:
            structured_data[wrcplkey][dTkey] = np.sqrt(np.sum([d_D["error"], d_N["error"]]))
            #structured_data[wrcplkey][dTkey] = d_D["error"]
        except TypeError:
            structured_data[wrcplkey][dTkey] = '-'

    for wr, cpl in itertools.product(evaluated_wr, coupling_schemes):
        wrcplkey = "WR{wr}, cpl={cpl}".format(wr=wr, cpl=cpl)
        wrkey = "WR{wr}".format(wr=wr)
        structured_data[wrcplkey]
        table.append([wrkey] + [structured_data[wrcplkey][dtkey] for dtkey in keys])

    return tabulate(table[1:], headers = table[0],tablefmt="latex_booktabs", floatfmt=".2e"), structured_data, keys


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Postprocessing of results')
    parser.add_argument('--prefix', '-p', help='Path prefix for results', type=str, default=os.path.join(os.path.dirname(os.path.realpath(__file__)),'../experiments'))
    args = parser.parse_args()

    #evaluated_wr = [11, 12, 13, 15,      21, 22, 23, 25, 31, 32, 33, 35, 51, 52, 53, 55, 101, 102, 103, 105, 110, 210, 310, 510, 1010]
    evaluated_wr = [11, 12,      15, 21, 22,     25,            51, 52,     55                                   ]
    evaluated_dT = [1.0, 0.5, 0.2, 0.1]
    #evaluated_dT = [1.0, 0.5]
    #coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]
    coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

    prefix = args.prefix

    table, data, keys = create_error_table(prefix, evaluated_wr, evaluated_dT, coupling_schemes)

    print(table)

    import itertools
    marker = itertools.cycle((',', '+', '.', 'o', '*'))

    for wr, cpl in itertools.product(evaluated_wr, coupling_schemes):
        wrcplkey = "WR{wr}, cpl={cpl}".format(wr=wr, cpl=cpl)
        plt.loglog(evaluated_dT, [data[wrcplkey][dtkey] for dtkey in keys], label=wrcplkey, marker=next(marker))

    plt.legend()
    plt.grid()
    plt.show()
