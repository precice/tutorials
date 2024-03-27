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
split_foldername = filename.split('/')
casename = split_foldername[0]
split_casename = casename.split('-')
solver = split_casename[-1]

if solver == 'python':
    df = pd.read_csv(filename, delimiter=';')
    if args.plotType == PlotType.U_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['position'].to_numpy())
        plt.title(PlotType.U_OVER_T.value)
    elif args.plotType == PlotType.V_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['velocity'].to_numpy())
        plt.title(PlotType.V_OVER_T.value)
    elif args.plotType == PlotType.TRAJECTORY.name:
        plt.plot(df['position'].to_numpy(), df['velocity'].to_numpy())
        plt.scatter([df['position'][0]], [df['velocity'][0]], label=f"(u,v) at t={df['time'][0]}")
        plt.scatter([df['position'].iloc[-1]], [df['velocity'].iloc[-1]],
                    label=f"(u,v) at t={df['time'].iloc[-1]}", marker="*")
        plt.title(PlotType.TRAJECTORY.value)
        plt.legend()
elif solver == 'fmi':
    df = pd.read_csv(filename, delimiter=',')
    if args.plotType == PlotType.U_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['mass.u'].to_numpy())
        plt.title(PlotType.U_OVER_T.value)
    elif args.plotType == PlotType.V_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['mass.v'].to_numpy())
        plt.title(PlotType.V_OVER_T.value)
    elif args.plotType == PlotType.TRAJECTORY.name:
        plt.plot(df[('mass.u')].to_numpy(), df['mass.v'].to_numpy())
        plt.scatter([df['mass.u'][0]], [df['mass.v'][0]], label=f"(u,v) at t={df['time'][0]}")
        plt.scatter([df['mass.u'].iloc[-1]], [df['mass.v'].iloc[-1]],
                    label=f"(u,v) at t={df['time'].iloc[-1]}", marker="*")
        plt.title(PlotType.TRAJECTORY.value)
        plt.legend()

plt.show()
