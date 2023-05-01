#! /usr/bin/env python3


from nutils import function, mesh, cli, solver, export
import treelog as log
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running utils")

    timestamp = 0.0

    # define the Nutils mesh
    nx = 60
    ny = 20

    grid = np.linspace(0, 6, nx + 1), np.linspace(0, 2, ny + 1)
    domain, geom = mesh.rectilinear(grid)
    domain = domain.withboundary(inflow="left", outflow="right", wall="top,bottom")

    # cloud of Gauss points
    gauss = domain.sample("uniform", degree=1)

    # Nutils namespace
    ns = function.Namespace(fallback_length=2)
    ns.x = geom

    # preCICE setup
    interface = precice.Interface("Heat", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Temperature-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    vertices = gauss.eval(ns.x)
    #print(vertices)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    # define "normal" state and "heating" state
    # values have 1 value per point, so it should have the same dim as the vertices' "x" dimension
    values_A = np.full(vertices.shape[0], 350.0)
    values_B = np.full(vertices.shape[0], 500.0)
    # 4 points around the middle
    #for i in range(640, 740):
    #    values_A[i] = 400.0
    #values_B[2440] = 400.0
    #values_B[2441] = 400.0
    #values_B[2442] = 400.0
    #values_B[2443] = 400.0

    # coupling data
    
    temperature_id = interface.get_data_id("Temperature", mesh_id)
    #temperature_id = 0

    precice_dt = interface.initialize()

    timestep = 0
    dt = 0.005

    # initialize the velocity values
    #temperature_values = np.zeros_like(vertices)
    #if timestep > 0.3 and timestep < 0.7:
    #    temperature_values = values_B
    #else:
    #    temperature_values = values_A

    
    while interface.is_coupling_ongoing():

        # TODO fix visualization
        #if timestep % 1 == 0:  # visualize
        #    print("in vis")
        #    bezier = domain.sample("bezier", 2)
        #    x, u = bezier.eval(["x_i"]@ ns), temperature_values
        #    print("x and u")
        #    with log.add(log.DataLog()):
        #        print("before")
        #        export.vtk("Heat_" + str(timestep), bezier.tri, x, T=u)
        #        print("after")

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)
        
        # if timestep between 0.3 and 0.7 write verB (heated), otherwise write verA (default)
        if interface.is_write_data_required(dt):
            if timestamp > 0.3 and timestamp < 2.3:
                interface.write_block_scalar_data(temperature_id, vertex_ids, values_B)
            else:
                interface.write_block_scalar_data(temperature_id, vertex_ids, values_A)
        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        timestamp += dt

    interface.finalize()


if __name__ == "__main__":
    cli.run(main)
