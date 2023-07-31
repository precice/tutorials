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
    uniform_grid = domain.sample("uniform", degree=1)

    # Nutils namespace
    ns = function.Namespace(fallback_length=2)
    ns.x = geom

    # preCICE setup
    interface = precice.Interface("Dummy-Velocity", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Dummy-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    vertices = uniform_grid.eval(ns.x)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    #values_A = np.full(vertices.shape[0], 350.0)
    values_B = np.full((vertices.shape[0], 2), [11.0, 0.0])

    # coupling data
    field_id = interface.get_data_id("Velocity", mesh_id)

    precice_dt = interface.initialize()

    timestep = 0
    dt = 0.005
    
    while interface.is_coupling_ongoing():

        # TODO fix visualization
        # if timestep % 1 == 0:  # visualize
        #     bezier = domain.sample("bezier", 2)
        #     x, u = bezier.eval(["x_i"]@ ns)[0], values_B
        #     with log.add(log.DataLog()):
        #         export.vtk("Dummy-Velocity_" + str(timestep), bezier.tri, x, T=u)

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)
        
        if interface.is_write_data_required(dt):
            #if timestamp > 0.3 and timestamp < 2.3:
            #    interface.write_block_scalar_data(field_id, vertex_ids, values_B)
            #else:
                interface.write_block_vector_data(field_id, vertex_ids, values_B)
        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        timestamp += dt

    interface.finalize()


if __name__ == "__main__":
    cli.run(main)
