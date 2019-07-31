#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:41:38 2019

@author: richyrich
"""

import numpy as np
import precice

configuration_file_name = "precice-config-dummy.xml"
participant_name = "Fluid"
mesh_name = "Fluid-Mesh-Nodes"

n = 1

solver_process_index = 0
solver_process_size = 1

interface = precice.Interface(participant_name, solver_process_index, solver_process_size)
interface.configure(configuration_file_name)
    
mesh_id = interface.get_mesh_id(mesh_name)
forces_id = interface.get_data_id("Forces0", mesh_id)

dimensions = interface.get_dimensions()
vertex = np.zeros(dimensions)
data_indices = np.zeros(n)

interface.set_mesh_vertices(mesh_id, n, vertex, data_indices)

dt = interface.initialize()
    
while interface.is_coupling_ongoing():
   
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        print("DUMMY: Writing iteration checkpoint")
        interface.fulfilled_action(precice.action_write_iteration_checkpoint())
    

    interface.write_block_vector_data(forces_id, n, data_indices, np.array((0,0,-1))) # write -1 for y-forces, 0 for x-forces
    dt = interface.advance(dt)
    
    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        print("DUMMY: Reading iteration checkpoint")
        interface.fulfilled_action(precice.action_read_iteration_checkpoint())
    else:
        print("DUMMY: Advancing in time")
    
interface.finalize()
print("DUMMY: Closing python solver dummy...")
