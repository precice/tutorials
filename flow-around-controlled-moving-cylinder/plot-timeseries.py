import pandas as pd
from matplotlib import pyplot as plt
import argparse
from enum import Enum


class PlotType(Enum):
    E_OVER_T = "error over time"
    P_OVER_T = "proportional term over time"
    I_OVER_T = "integral term over time"
    D_OVER_T = "derivative term over time"
    U1_OVER_T = "control output 1 over time"
    U2_OVER_T = "control output 2 over time"
    Y1_OVER_T = "control input 1 over time"
    Y2_OVER_T = "control input 2 over time"


parser = argparse.ArgumentParser()
parser.add_argument("csvFile", help="Path to CSV file.", type=str)
parser.add_argument(
    "plotType",
    help="Which data should be plotted, e.g. E_OVER_T = Error over time",
    type=str,
    choices=[
        pt.name for pt in PlotType])
args = parser.parse_args()

filename = args.csvFile
split_filename = filename.split('/')
solver = split_filename[0]
if solver == ".":
    solver = split_filename[1]


if solver == 'controller-fmi':
    df = pd.read_csv(filename, delimiter=',')
    if args.plotType == PlotType.E_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['e'].to_numpy())
        plt.title(PlotType.E_OVER_T.value)
    elif args.plotType == PlotType.P_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['P'].to_numpy())
        plt.title(PlotType.P_OVER_T.value)
    elif args.plotType == PlotType.I_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['I'].to_numpy())
        plt.title(PlotType.I_OVER_T.value)
    elif args.plotType == PlotType.D_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['D'].to_numpy())
        plt.title(PlotType.D_OVER_T.value)
    elif args.plotType == PlotType.U1_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['u_1'].to_numpy())
        plt.title(PlotType.U1_OVER_T.value)
    elif args.plotType == PlotType.U2_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['u_2'].to_numpy())
        plt.title(PlotType.U2_OVER_T.value)
    elif args.plotType == PlotType.Y1_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['y_1'].to_numpy())
        plt.title(PlotType.Y1_OVER_T.value)
    elif args.plotType == PlotType.Y2_OVER_T.name:
        plt.plot(df['time'].to_numpy(), df['y_2'].to_numpy())
        plt.title(PlotType.Y2_OVER_T.value)
    else:
        print("Warning: Controller can not plot this value. Please type run python3 plot-timeseries.py -h to see the available commands.")

plt.show()
