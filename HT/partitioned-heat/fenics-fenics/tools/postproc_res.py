import os
from tools.coupling_schemes import CouplingScheme
import itertools
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Postprocessing of results')
parser.add_argument('--evaluated_wr_left', '-wrl', type=int, default=1)
parser.add_argument('--evaluated_wr_right', '-wrr', type=int, default=1)
parser.add_argument('--evaluated_dT', '-dT', type=float, default=1.0)
parser.add_argument('--gamma', '-g', type=float, default=0.0)
parser.add_argument('--coupling_scheme', '-cpl', type=str, default='d')
parser.add_argument('--prefix', '-p', help='Path prefix for results', type=str, default='.')
args = parser.parse_args()

#evaluated_wr = [11, 12, 13, 15,      21, 22, 23, 25, 31, 32, 33, 35, 51, 52, 53, 55, 101, 102, 103, 105, 110, 210, 310, 510, 1010]
#evaluated_wr = [11, 12,      15, 110, 21, 22,     25, 210,            51, 52,     55,                                    510, 101, 102, 105, 1010]
#evaluated_dT = [1.0, 0.5, 0.2, 0.1]
#coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name, CouplingScheme.SERIAL_FIRST_NEUMANN.name, CouplingScheme.PARALLEL.name]

evaluated_wr_left = [args.evaluated_wr_left]
evaluated_wr_right = [args.evaluated_wr_right]
evaluated_dT = [args.evaluated_dT]
if args.coupling_scheme == 'd':
    coupling_schemes = [CouplingScheme.SERIAL_FIRST_DIRICHLET.name]
if args.coupling_scheme == 'n':
    coupling_schemes = [CouplingScheme.SERIAL_FIRST_NEUMANN.name]
if args.coupling_scheme == 'p':
    coupling_schemes = [CouplingScheme.PARALLEL.name]

prefix = args.prefix
"""
if args.gamma == 0.0:
    prefix = "/home/benjamin/2019paperqnwr/PartitionedHeatEquation/experiments_gamma_0_2019_7_13_afternoon/"
elif args.gamma == 1.0:
    prefix = "/home/benjamin/2019paperqnwr/PartitionedHeatEquation/experiments_gamma_1_2019_7_14_morning/"
"""

for wr_left, wr_right, dT, coupling_scheme in itertools.product(evaluated_wr_left, evaluated_wr_right, evaluated_dT, coupling_schemes):
    wr_tag = "WR{wr_left}{wr_right}".format(wr_left=wr_left, wr_right=wr_right)
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
    
    i = 0
    legend_handles = []
    legend_lables = []
    if args.coupling_scheme == 'd':
        first_quantity = 'q'
        second_quantity = 'T'
    if args.coupling_scheme == 'n':
        first_quantity = 'q'
        second_quantity = 'T'    

    for wr_it in range(wr_left):
        plt.semilogy(data_dict['total_it'],data_dict['resNorm({i})'.format(i=i)], label="Res: {first_quantity}{wr_it}".format(first_quantity = first_quantity, wr_it = wr_it+1))
        i+=1
    for wr_it in range(wr_right):
        plt.semilogy(data_dict['total_it'],data_dict['resNorm({i})'.format(i=i)], label="Res: {second_quantity}{wr_it}".format(second_quantity = second_quantity, wr_it = wr_it+1))

    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])    
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.title("{wr_tag}, {window_tag}, {cpl_tag}, gamma={gamma}".format(wr_tag=wr_tag, window_tag=window_tag, cpl_tag=coupling_scheme, gamma=args.gamma))
    fig_path = os.path.join(folder,convergence_file_name[:-4]+"_fig.png")
    plt.savefig(fig_path)
    print("output written to "+ fig_path)
