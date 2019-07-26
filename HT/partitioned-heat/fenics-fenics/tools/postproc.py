import os
from coupling_schemes import CouplingScheme
import itertools
from tabulate import tabulate
import argparse

parser = argparse.ArgumentParser(description='Postprocessing of results')
parser.add_argument('--prefix', '-p', help='Path prefix for results', type=str, default='..')
args = parser.parse_args()

#evaluated_wr = [11, 12, 13, 15,      21, 22, 23, 25, 31, 32, 33, 35, 51, 52, 53, 55, 101, 102, 103, 105, 110, 210, 310, 510, 1010]
evaluated_wr = [11, 12,      15, 110, 21, 22,     25, 210,            51, 52,     55,                                    510, 101, 102, 105, 1010]
evaluated_dT = [1.0, 0.5, 0.2, 0.1]
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]

simulationTime = 10.0
data = []
data.append(["WR", "dT", "cpl", "#steps", "#it", "#it/#steps"])
prefix = args.prefix

for wr, dT, coupling_scheme in itertools.product(evaluated_wr, evaluated_dT, coupling_schemes):
    wr_tag = "WR{wr}".format(wr=wr)
    window_tag = "dT{dT}".format(dT=dT)
    folder = os.path.join(prefix, "experiments", wr_tag, window_tag, coupling_scheme)
    try:
        with open(os.path.join(folder, "precice-HeatDirichlet-iterations.log")) as file:
            for line in file.readlines():
                lastline = line
            try:
                lastline = lastline.split("  ")
                number_of_windows = float(lastline[0])
                total_iterations = float(lastline[1])
                avg_its = total_iterations / number_of_windows
                if number_of_windows == simulationTime / dT:
                    data.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "#steps": number_of_windows, "#it": total_iterations, "#it/#steps": avg_its})
                else:
                    print("File {name} erroneous.".format(
                        name=os.path.join(folder, "precice-HeatDirichlet-iterations.log")))
                    data.append(
                        {"WR": wr, "dT": dT, "cpl": coupling_scheme, "#steps": "-", "#it": "-", "#it/#steps": "x"})
            except AttributeError:
                print("File {name} erroneous.".format(name=os.path.join(folder, "precice-HeatDirichlet-iterations.log")))
                data.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "#steps": "-", "#it": "-", "#it/#steps": "x"})
    except FileNotFoundError:
        data.append({"WR": wr, "dT": dT, "cpl": coupling_scheme, "#steps": "-", "#it": "-", "#it/#steps": "-"})

table = []
keys = ["dT{dT}".format(dT=dT) for dT in evaluated_dT]
table.append(["#it/#steps"] + keys)

structured_data = dict()

for d in data[1:]:
    wrcplkey = "WR{wr}, cpl={cpl}".format(wr=d["WR"], cpl=d["cpl"])
    try:
        structured_data[wrcplkey]
    except KeyError:
        structured_data[wrcplkey] = dict()
    dTkey = "dT{dT}".format(dT=d["dT"])
    structured_data[wrcplkey][dTkey] = d["#it/#steps"]
    
for wr, cpl in itertools.product(evaluated_wr, coupling_schemes):
    wrcplkey = "WR{wr}, cpl={cpl}".format(wr=wr, cpl=cpl)
    structured_data[wrcplkey]
    table.append([wrcplkey] + [structured_data[wrcplkey][dtkey] for dtkey in keys])

print(tabulate(table[1:], headers = table[0],tablefmt="pipe"))
