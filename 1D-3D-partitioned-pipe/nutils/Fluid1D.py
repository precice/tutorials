#!/usr/bin/env python3

from numpy.core.fromnumeric import reshape
from evtk.hl import pointsToVTK
import numpy as np
import precice

def main():
    
    # number of nodes, length of domain and space interval

    length = 10
    n = 20
    h = 1 / n
    t = 0
    counter = 0

    # generate mesh

    x = np.zeros(n+1)
    y = np.zeros(n+1)
    z = np.linspace(0,length,n+1)

    x_interface = np.linspace(-5,5,20)
    y_interface = np.linspace(-5,5,20)

    # initial data

    u = np.zeros((n+1,3))
    p = np.zeros(n+1)
    rhs = np.zeros(n+1)

    # preCICE setup

    participant_name = "Fluid1"
    config_file_name = "./precice-config.xml"
    solver_process_index = 0
    solver_process_size = 1
    interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)

    mesh_name = "Fluid1-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)

    velocity_name = "Velocity"
    velocity_id = interface.get_data_id(velocity_name, mesh_id)

    pressure_name = "Pressure"
    pressure_id = interface.get_data_id(pressure_name, mesh_id)

    positions = [0, 0, z[-1]]

    vertex_ids = interface.set_mesh_vertex(mesh_id, positions)

    precice_dt = interface.initialize()


    while interface.is_coupling_ongoing():
        if interface.is_action_required(
            precice.action_write_iteration_checkpoint()):

            u_iter = u
            p_iter = p
            t_iter = t
            counter_iter = counter

            interface.mark_action_fulfilled(
            precice.action_write_iteration_checkpoint())

        # determine time step size

        dt = 0.025
        dt = np.minimum(dt,precice_dt)

        # set boundary conditions

        u[0,2] = 1
        # u[1,2] = 1  # dirichlet velocity inlet

        u[-1,2] = u[-2,2]   # neumann velocity outlet

        p[0] = p[1] # neumann pressure inlet
        # p[-1] = 0   # dirichlet pressure outlet
        if interface.is_read_data_available():  # get dirichlet pressure outlet value from 3D solver
            p_read_in = interface.read_scalar_data(pressure_id, vertex_ids)
            p[-1] = p_read_in

        # compute right-hand side of 1D PPE

        for i in range(n-1):

            rhs[i+1] = (1 / dt) * ((u[i+2,2] - u[i,2]) / 2*h)

        # solve the PPE using a SOR solver

        tolerance = 0.001
        error = 1
        omega = 1.7
        max_iter = 1000
        iter = 0

        while error >= tolerance:

            p[0] = p[1] # renew neumann pressure inlet
            
            for i in range(n-1):
                p[i+1] = (1-omega) * p[i+1] + ((omega * h**2) / 2) * (((p[i] + p[i+2]) / h**2) - rhs[i+1])

            sum = 0
            for i in range(n-1):
                val = ((p[i] - 2*p[i+1] + p[i+2]) / h**2) - rhs[i+1]
                sum += val*val

            error = np.sqrt(sum/n)

            iter += 1
            if iter >= max_iter:
                print("SOR solver did not converge.\n")
                break
                

        # calculate new velocities

        for i in range(n-1):
            u[i+1,2] = u[i+1,2] - dt * ((p[i+2] - p[i+1]) / h)


        # write velocity to 3D solver

        write_vel = u[-2,:]

        if interface.is_write_data_required(dt):    # write new velocities to 3D solver
            interface.write_vector_data(
            velocity_id, vertex_ids, write_vel)

        # transform data and write output data to vtk files

        u_print = np.reshape(u[:,2], n+1)
        p_print = np.reshape(p, n+1)

        u_print = np.ascontiguousarray(u_print, dtype=np.float32)
        p_print = np.ascontiguousarray(p_print, dtype=np.float32)
        
        filename = "./results/1D/Fluid1D_" + str(counter)
        pointsToVTK(filename, x, y, z, data = {"U" : u_print, "p" : p_print})

        # advance simulation time

        dt = interface.advance(dt)
        t = t + dt
        counter += 1

        if interface.is_action_required(
            precice.action_read_iteration_checkpoint()):

            u = u_iter
            p = p_iter
            t = t_iter
            counter = counter_iter

            interface.mark_action_fulfilled(
            precice.action_read_iteration_checkpoint())

    interface.finalize()


if __name__ == "__main__":
    main()