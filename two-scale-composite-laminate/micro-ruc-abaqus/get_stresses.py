"""
This file is used to extract stresses from the odb file of the RUC simulation and write them to a file called stresses.txt.
The file stresses.txt is then read by the ruc_abaqus_restart.py file to get the stresses for the current iteration.
File written by Minh Hoang Nguyen, mhoangn@umich.edu, 2023
"""

from abaqus import *
from abaqusConstants import *
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import numpy as np
import os

import __main__

working_dir = os.getcwd()

length = len(working_dir)
id_as_string = working_dir[length - 2:]

ruc_full_path = working_dir + '/RUC_' + id_as_string + '.odb'

# Check that the odb file exists
assert os.path.exists(ruc_full_path), "RUC " + id_as_string + " does not have an .odb file at location: " + working_dir

o1 = session.openOdb(name=ruc_full_path, readOnly=False)
odb = session.odbs[ruc_full_path]

labels = ['X', 'Y', 'Z']

for iRefP in range(3 + 1)[1:]:
    label = labels[iRefP - 1]
    for iDir in range(3 + 1)[1:]:
        data_name = 'RF{iDir} PI: rootAssembly N: {iRefP} NSET SETRP{label}-1'.format(
            iDir=iDir, iRefP=iRefP, label=label)
        outputVariable_name = 'Reaction force: RF{iDir} PI: rootAssembly Node {iRefP} in NSET SETRP{label}'.format(
            iDir=iDir, iRefP=iRefP, label=label)
        session.XYDataFromHistory(odb=odb, name=data_name, outputVariableName=outputVariable_name,
                                  steps=('Step-1', ), __linkedVpName__='Viewport: 1')

# Areas
stresses = np.zeros(6)  # Only the following components of the stress tensor are non-zero: [11, 22, 33, 12, 23, 13]

Lx = 2.600E-03
Ly = 6.447E-03
Lz = 11.167E-03

A1 = Ly * Lz
A2 = Lx * Lz
A3 = Lx * Ly

x = np.array(session.xyDataObjects['RF1 PI: rootAssembly N: 1 NSET SETRPX-1'])
stresses[0] = x[1][1] / A1  # Sig11

x = np.array(session.xyDataObjects['RF2 PI: rootAssembly N: 2 NSET SETRPY-1'])
stresses[1] = x[1][1] / A2  # Sig22

x = np.array(session.xyDataObjects['RF3 PI: rootAssembly N: 3 NSET SETRPZ-1'])
stresses[2] = x[1][1] / A3  # Sig33

# Tau12
x = np.array(session.xyDataObjects['RF1 PI: rootAssembly N: 2 NSET SETRPY-1'])
stresses[3] = x[1][1] / A2  # Tau12

# Tau23
x = np.array(session.xyDataObjects['RF3 PI: rootAssembly N: 2 NSET SETRPY-1'])
stresses[4] = x[1][1] / A2  # Tau23

# Tau13
x = np.array(session.xyDataObjects['RF1 PI: rootAssembly N: 3 NSET SETRPZ-1'])
stresses[5] = x[1][1] / A3  # Tau13

# Convert stresses to strings
stresses_as_strings = []
for i in range(6):
    stresses_as_strings.append('{}\n'.format(stresses[i]))

# Create a file and write stresses to it
output_file = open('stresses.txt', 'w')

output_file.writelines(stresses_as_strings)
output_file.close()
