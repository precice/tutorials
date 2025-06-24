import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

files = {"default": "chemistry-fenics/output/chemical_out.csv"}
datasets = {}

headers = ['Time', 'A', 'B', 'C']

for key, path in files.items():
    datasets[key] = pd.read_csv(files[key], names=headers, delimiter=' ')
    plt.plot(datasets[key]['Time'], datasets[key]['A'], label=key)
    plt.title("Total amount of A")
    plt.xlabel("Time")
    plt.ylabel("Total amount")
    plt.legend()
    plt.show()

for key, path in files.items():
    datasets[key] = pd.read_csv(files[key], names=headers, delimiter=' ')
    plt.plot(datasets[key]['Time'], datasets[key]['B'], label=key)
    plt.title("Total amount of B")
    plt.xlabel("Time")
    plt.ylabel("Total amount")
    plt.legend()
    plt.show()

for key, path in files.items():
    datasets[key] = pd.read_csv(files[key], names=headers, delimiter=' ')
    plt.plot(datasets[key]['Time'], datasets[key]['C'], label=key)
    plt.title("Total amount of C")
    plt.xlabel("Time")
    plt.ylabel("Total amount")
    plt.legend()
    plt.show()
