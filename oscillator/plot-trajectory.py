import pandas as pd
from matplotlib import pyplot as plt
import argparse
from enum import Enum


class PlotType(Enum):
    U_OVER_T = "position over time"
    V_OVER_T = "velocity over time"
    TRAJECTORY = "velocity over position (trajectory)"


parser = argparse.ArgumentParser()
parser.add_argument("csvFile", help="CSV file.", type=str)
parser.add_argument("plotType", help="Plot type.", type=str, choices=[pt.name for pt in PlotType])
args = parser.parse_args()

filename = args.csvFile
split_filename = filename.split('/')
solver = split_filename[0]

if solver == 'python':
    df = pd.read_csv(filename, delimiter=';')
    if args.plotType == PlotType.U_OVER_T.name:
        plt.plot(df['time'], df['position'])
        plt.title(PlotType.U_OVER_T.value)
    elif args.plotType == PlotType.V_OVER_T.name:
        plt.plot(df['time'], df['velocity'])
        plt.title(PlotType.V_OVER_T.value)
    elif args.plotType == PlotType.TRAJECTORY.name:
        plt.plot(df['position'], df['velocity'])
        plt.scatter([df['position'][0]], [df['velocity'][0]], label=f"(u,v) at t={df['time'][0]}")
        plt.scatter([df['position'].iloc[-1]], [df['velocity'].iloc[-1]],
                    label=f"(u,v) at t={df['time'].iloc[-1]}", marker="*")
        plt.title(PlotType.TRAJECTORY.value)
        plt.legend()   
elif solver == 'fmi':
    df = pd.read_csv(filename, delimiter=',')
    if args.plotType == PlotType.U_OVER_T.name:
        plt.plot(df['time'], df['mass.u'])
        plt.title(PlotType.U_OVER_T.value)
    elif args.plotType == PlotType.V_OVER_T.name:
        plt.plot(df['time'], df['mass.v'])
        plt.title(PlotType.V_OVER_T.value)
    elif args.plotType == PlotType.TRAJECTORY.name:
        plt.plot(df[('mass.u')], df['mass.v'])
        plt.scatter([df['mass.u'][0]], [df['mass.v'][0]], label=f"(u,v) at t={df['time'][0]}")
        plt.scatter([df['mass.u'].iloc[-1]], [df['mass.v'].iloc[-1]],
                    label=f"(u,v) at t={df['time'].iloc[-1]}", marker="*")
        plt.title(PlotType.TRAJECTORY.value)
        plt.legend()

plt.show()
