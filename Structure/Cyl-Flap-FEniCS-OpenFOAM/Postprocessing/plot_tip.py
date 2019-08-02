#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 15:07:54 2019

@author: richyrich
"""

import numpy as np
import matplotlib.pyplot as plt

fenics_watchpoint = np.loadtxt('../precice-fenics-watchpoint-point1.log', skiprows = 1)
ccx_watchpoint = np.loadtxt('/home/richyrich/Bachelorarbeit/OpenFOAM-CalculiX/precice-Fluid-watchpoint-point1.log', skiprows = 1)


plt.plot(fenics_watchpoint[:,0],fenics_watchpoint[:,9])
