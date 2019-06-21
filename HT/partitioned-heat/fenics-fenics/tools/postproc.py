import os
from participants import Participant
import itertools

evaluated_wr = [11, 12, 13, 15, 21, 22, 23, 25, 31, 32, 33, 35, 51, 52, 53, 55, 101, 102, 103, 105, 110, 210, 310, 510, 1010]
evaluated_dT = [1.0, 0.5, 0.2, 0.1]
first_participants = [Participant.DIRICHLET.name, Participant.NEUMANN.name]

data = []
data.append(["WR", "dT", "first", "#steps", "#it", "#it/#steps"])

for wr, dT, first_participant in itertools.product(evaluated_wr, evaluated_dT, first_participants):
    print(wr)
    print(dT)
    print(first_participant)
    wr_tag = "WR{wr}".format(wr=wr)
    window_tag = "dT{dT}".format(dT=dT)
    first_participant_tag = "first_{}".format(first_participant)
    folder = os.path.join("../experiments", wr_tag, window_tag, first_participant_tag)
    
    try:
        with open(os.path.join(folder, "precice-HeatDirichlet-iterations.log")) as file:
            for line in file.readlines():
                lastline = line
            lastline = lastline.split("  ")
            number_of_windows = float(lastline[0])
            total_iterations = float(lastline[1])
            avg_its = total_iterations / number_of_windows
            data.append({"WR": wr, "dT": dT, "first": first_participant, "#steps": number_of_windows, "#it": total_iterations, "#it/#steps": avg_its})
    except FileNotFoundError:
        data.append({"WR": wr, "dT": dT, "first": first_participant, "#steps": "-", "#it": "-", "#it/#steps": "-"})

table = []
keys = ["dT{dT}".format(dT=dT) for dT in evaluated_dT]
table.append(["#it/#steps"] + keys)

structured_data = dict()

for d in data[1:]:
    print(d)
    wrfirstkey = "WR{wr}, first={first}".format(wr=d["WR"], first=d["first"])
    print(wrfirstkey)
    try:
        structured_data[wrfirstkey]
    except KeyError:
        structured_data[wrfirstkey] = dict()
    dTkey = "dT{dT}".format(dT=d["dT"])
    structured_data[wrfirstkey][dTkey] = d["#it/#steps"]
    
for wr, first in itertools.product(evaluated_wr, first_participants):
    wrfirstkey = "WR{wr}, first={first}".format(wr=wr, first=first)
    structured_data[wrfirstkey]
    table.append([wrfirstkey] + [structured_data[wrfirstkey][dtkey] for dtkey in keys])

for line in table:
    print(line)

from tabulate import tabulate

print(tabulate(table[1:], headers = table[0]))
