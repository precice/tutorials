#!/usr/bin/python

import vtk
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

T = 100  # number of timesteps performed

arrayname = sys.argv[1]  # Which dataset should be plotted?
data_path = sys.argv[2]  # Where is the data?


def file_name_generator(id): return data_path + str(id) + ".vtk"


print("reading data from array with name = %s" % arrayname)
print("parsing datasets named %s*.vtk" % data_path)

values_for_all_t = T * [None]


for t in range(T):

    # read the vtk file as an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file_name_generator(t))
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    # parse the data
    grid = reader.GetOutput()
    point_data = grid.GetPointData().GetArray(arrayname)
    points = grid.GetPoints()
    N = grid.GetNumberOfPoints()  # How many gridpoints do exist?

    if point_data is None:  # check if array exists in dataset
        print("array with name %s does not exist!" % arrayname)
        print("exiting.")
        quit()

    value_at_t = []
    spatial_mesh = []

    n = point_data.GetNumberOfComponents()

    for i in range(N):  # parse data from vtk array into list

        x, y, z = grid.GetPoint(i)  # read coordinates of point
        spatial_mesh += [x]  # only store x component

        v = np.zeros(n)  # initialize empty butter array
        point_data.GetTuple(i, v)  # read value into v
        value_at_t += [np.linalg.norm(v)]

    values_for_all_t[t] = value_at_t


values_for_all_t = np.array(values_for_all_t)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(spatial_mesh, range(T))

# uncomment depending on what quantity you want to plot
ax.plot_surface(X, Y, values_for_all_t, cmap='viridis', edgecolor='black')
plt.xlabel("space")
plt.ylabel("time")
plt.title(arrayname + " from " + data_path + "*.vtk")
plt.show()
