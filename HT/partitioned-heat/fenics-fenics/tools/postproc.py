import os
from tools.coupling_schemes import CouplingScheme
import itertools
from tabulate import tabulate
import argparse


def create_qn_table(prefix, evaluated_wr, evaluated_dT, coupling_schemes, simulation_time=10):

    data = []
    data.append(["WR", "dT", "cpl", "#steps", "#it", "#it/#steps"])
    print(evaluated_dT)
    for wr, dT, coupling_scheme in itertools.product(evaluated_wr, evaluated_dT, coupling_schemes):
        wr_tag = "WR{wr}".format(wr=wr)
        window_tag = "dT{dT}".format(dT=dT)
        folder = os.path.join(prefix, wr_tag, window_tag, coupling_scheme)
        try:
            with open(os.path.join(folder, "precice-HeatDirichlet-iterations.log")) as file:
                for line in file.readlines():
                    lastline = line
                try:
                    lastline = lastline.split("  ")
                    number_of_windows = float(lastline[0])
                    total_iterations = float(lastline[1])
                    avg_its = total_iterations / number_of_windows
                    if number_of_windows == simulation_time / float(dT):
                        data.append({"WR": wr, "dT": dT, "#steps": number_of_windows, "#it": total_iterations, "#it/#steps": avg_its})
                    else:
                        print("File {name} erroneous.".format(
                            name=os.path.join(folder, "precice-HeatDirichlet-iterations.log")))
                        data.append(
                            {"WR": wr, "dT": dT, "#steps": "-", "#it": "-", "#it/#steps": "x"})
                except (AttributeError, ValueError):
                    print("File {name} erroneous.".format(name=os.path.join(folder, "precice-HeatDirichlet-iterations.log")))
                    data.append({"WR": wr, "dT": dT, "#steps": "-", "#it": "-", "#it/#steps": "x"})
        except FileNotFoundError:
            data.append({"WR": wr, "dT": dT, "#steps": "-", "#it": "-", "#it/#steps": "-"})

    table = []
    keys = [ "WR{wr}".format(wr=wr) for wr in evaluated_wr]
    table.append(["dT"] + keys)

    structured_data = dict()

    for d in data[1:]:
        wrkey = "WR{wr}".format(wr=d["WR"])
        try:
            structured_data[wrkey]
        except KeyError:
            structured_data[wrkey] = dict()
        dTkey = d["dT"]
        structured_data[wrkey][dTkey] = d["#it/#steps"]

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
    csv_file = "iterations_raw.csv"
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for data in data_dict:
            writer.writerow(data)

    return tabulate(table[1:], headers=table[0], tablefmt="latex_booktabs", floatfmt=".2f"), structured_data, keys


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

    table, data, keys = create_qn_table(prefix, evaluated_wr, evaluated_dT, coupling_schemes)

