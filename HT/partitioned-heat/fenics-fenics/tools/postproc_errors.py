import os
from tools.coupling_schemes import CouplingScheme
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
                    data_D.append({"WR": wr, "dT": dT, "error": error})
            except FileNotFoundError:
                data_D.append({"WR": wr, "dT": dT, "error": "-"})
        else:
            data_D.append({"WR": wr, "dT": dT, "error": "-"})

        assert(glob.glob(path_N).__len__() <= 1)
        if(glob.glob(path_N).__len__() == 1):
            try:
                with open(glob.glob(path_N)[0]) as file:
                    for line in file.readlines():
                        lastline = line
                    lastline = lastline.split(", ")
                    dt = float(lastline[0])
                    error = float(lastline[1])
                    data_N.append({"WR": wr, "dT": dT, "error": error})
            except FileNotFoundError:
                data_N.append({"WR": wr, "dT": dT, "error": "-"})
        else:
            data_N.append({"WR": wr, "dT": dT, "error": "-"})

    table = []
    keys = [ "WR{wr}".format(wr=wr) for wr in evaluated_wr]
    table.append(["dT"] + keys)

    structured_data = dict()

    for d_D, d_N in zip(data_D[1:],data_N[1:]):
        wrkey = "WR{wr}".format(wr=d_D["WR"])

        try:
            structured_data[wrkey]
        except KeyError:
            structured_data[wrkey] = dict()
        dTkey = d_D["dT"]
        try:
            structured_data[wrkey][dTkey] = float(np.sqrt(np.sum([d_D["error"], d_N["error"]])))
        except TypeError:
            structured_data[wrkey][dTkey] = np.nan

    dTkeys = [dT for dT in evaluated_dT]
    for dTkey in dTkeys:
        table.append([dTkey] + [structured_data[key][dTkey] for key in keys])

    data_dict = []       
    for dT in evaluated_dT:
        dict_line = dict()
        dict_line['dT'] = dT
        for key in structured_data.keys():
            dict_line[key] = structured_data[key][dT]
        data_dict.append(dict_line)

    import csv
    headers = ['dT'] + list(structured_data.keys())
    csv_file = "errors_raw.csv"
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for data in data_dict:
            writer.writerow(data)


    return tabulate(table[1:], headers = table[0],tablefmt="latex_booktabs", floatfmt=".2e"), structured_data, dTkeys


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Postprocessing of results')
    parser.add_argument('--prefix', '-p', help='Path prefix for results', type=str, default=os.path.join(os.path.dirname(os.path.realpath(__file__)),'../experiments'))
    args = parser.parse_args()

    #evaluated_wr = [11, 12, 13, 15,      21, 22, 23, 25, 31, 32, 33, 35, 51, 52, 53, 55, 101, 102, 103, 105, 110, 210, 310, 510, 1010]
    evaluated_wr = [11, 12,      15, 21, 22,     25,            51, 52,     55                                   ]
    evaluated_dT = [5.0, 2.0, 1.0, 0.5, 0.2]
    #evaluated_dT = [1.0, 0.5]
    #coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]
    coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

    prefix = args.prefix

    table, data, dTkeys = create_error_table(prefix, evaluated_wr, evaluated_dT, coupling_schemes)

    import itertools
    marker = itertools.cycle((',', '+', '.', 'o', '*'))

    for key in data.keys():
        print("{}:{}".format(key,data[key]))
        plt.loglog(evaluated_dT, [data[key][dTkey] for dTkey in dTkeys], label=key, marker=next(marker))

    plt.legend()
    plt.grid()
    plt.show()
