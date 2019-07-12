import os
from coupling_schemes import CouplingScheme
import itertools
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Postprocessing of results')
parser.add_argument('--prefix', '-p', help='Path prefix for results', type=str, default='..')
args = parser.parse_args()

#evaluated_wr = [11, 12, 13, 15,      21, 22, 23, 25, 31, 32, 33, 35, 51, 52, 53, 55, 101, 102, 103, 105, 110, 210, 310, 510, 1010]
#evaluated_wr = [11, 12,      15, 110, 21, 22,     25, 210,            51, 52,     55,                                    510, 101, 102, 105, 1010]
#evaluated_dT = [1.0, 0.5, 0.2, 0.1]
#coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]

evaluated_wr = [12]
evaluated_dT = [1.0]
coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]

prefix = args.prefix

for wr, dT, coupling_scheme in itertools.product(evaluated_wr, evaluated_dT, coupling_schemes):
    wr_tag = "WR{wr}".format(wr=wr)
    window_tag = "dT{dT}".format(dT=dT)
    
    if coupling_scheme == CouplingScheme.SERIAL_FIRST_DIRICHLET.name:
        convergence_file_name = "precice-HeatNeumann-convergence.log"
    elif coupling_scheme == CouplingScheme.SERIAL_FIRST_NEUMANN.name:
        convergence_file_name = "precice-HeatDirichlet-convergence.log"
    
    folder = os.path.join(prefix, "experiments", wr_tag, window_tag, coupling_scheme)
    with open(os.path.join(folder, convergence_file_name)) as file:            
        lines = file.readlines()
        headers = lines[0].split("  ")
        headers.append("total_it")
        data_dict = dict.fromkeys(headers)
        for k in data_dict.keys():
            data_dict[k] = list()
        total_iteration = 0
        for line in lines[1::]:
            total_iteration += 1
            i = 0
            print(line.split("  ")[:-1])
            for item in line.split("  ")[:-1]:
                print("{key}:{value}".format(key=headers[i], value=float(item)))
                try:
                    data_dict[headers[i]].append(float(item))
                except ValueError:
                    print(item)
                    data_dict[headers[i]].append("INVALID")
                i+=1
            print("{key}:{value}".format(key="total_it", value=total_iteration))
            data_dict["total_it"].append(total_iteration)               
            print(data_dict)
        print(headers)
    print(data_dict[headers[1]])

    plt.semilogy(data_dict['total_it'],data_dict['resNorm(0)'])
    plt.semilogy(data_dict['total_it'],data_dict['resNorm(1)'])
    plt.semilogy(data_dict['total_it'],data_dict['resNorm(2)'])
    #plt.plot(data_dict['total_it'],data_dict['Iteration'])
    plt.show()

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
