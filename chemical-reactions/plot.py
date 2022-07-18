import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

files = {"nn32" : 'out_nn32.csv',  "nn64" : 'out_nn64.csv', "lci32" : 'out_lci32.csv', "lci50" : 'out_lci50.csv'}
datasets = {}

headers = ['Time', 'A', 'B', 'C']

for key, path in files.items():
    datasets[key] = pd.read_csv(files[key], names=headers, delimiter = ' ')
    plt.plot(datasets[key]['Time'], datasets[key]['C'], label = key)

plt.legend()
plt.show()