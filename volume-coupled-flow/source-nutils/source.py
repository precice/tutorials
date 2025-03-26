#! /usr/bin/env python3


from nutils import function, mesh, cli, solver, export
import treelog as log
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running Nutils")

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
    participant = precice.Participant("Source-Velocity", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Source-Mesh"
    vertices = uniform_grid.eval(ns.x)
    vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

    source_values = np.full((vertices.shape[0], 2), [10.0, 0.0])

    # coupling data
    data_name = "Velocity"

    participant.initialize()

    timestep = 0
    solver_dt = 0.005

    while participant.is_coupling_ongoing():

        precice_dt = participant.get_max_time_step_size()

        # potentially adjust non-matching timestep sizes
        dt = min(solver_dt, precice_dt)

        participant.write_data(mesh_name, data_name, vertex_ids, source_values)
        # do the coupling
        participant.advance(dt)

        # advance variables
        timestep += 1
        timestamp += dt

    participant.finalize()


if __name__ == "__main__":
    cli.run(main)
