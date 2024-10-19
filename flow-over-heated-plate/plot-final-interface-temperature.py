#!/usr/bin/env python3
import vtk
from matplotlib import pyplot as plt
import numpy as np
import os

def vtk_to_dict(case):
    vtkFileName = "solid-{}/precice-exports/Fluid-Mesh-Solid.dt100.vtk".format(case)
    if not os.path.exists(vtkFileName):
        print("No file found for " + vtkFileName)
        return {} # return empty dict if file not found
    
    # read the vtk file as an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFileName)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    # obtain the data
    data = reader.GetOutput()
    n_data = data.GetPointData().GetNumberOfTuples()

    data_dict = {}

    for i in range(n_data):
        data_dict[data.GetPoint(i)] = data.GetPointData().GetArray("Temperature").GetValue(i)
    return data_dict



def main():
    case_labels = {
        'fenics': 'Fluid-FEniCS',
        'openfoam': 'Fluid-OpenFOAM',
        'nutils': 'Fluid-Nutils',
        'dunefem': 'Fluid-DuneFem'}
    styles = [':', '-', '--']
    colors = ['r', 'b', 'g', 'k']

    for i, case in enumerate(case_labels.keys()):
        case_data = vtk_to_dict(case)
        if not case_data:
            continue
        x, t = [p[0] for p in case_data.keys()], np.array(list(case_data.values()))

        # sort by x
        combined = list(zip(x,t))
        combined.sort()
        x, t = zip(*combined)
        x = np.array(x)
        t = np.array(t)

        theta = (t - 300) / (310 - 300)
        plt.plot(x, theta, colors[i % 4] + styles[i % 3], label=case_labels[case])

    plt.ylabel("Theta")
    plt.xlabel("x-coordinate along coupling interface")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()