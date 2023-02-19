#! /usr/bin/env python3

#
# Advection-Diffusion equation for a single species with a velocity field read from preCICE on the complete volume.
#

from nutils import function, mesh, cli, solver, export
import treelog as log
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running utils")

    # define the Nutils mesh
    nx = 60
    ny = 20

    grid = np.linspace(0, 6, nx + 1), np.linspace(0, 2, ny + 1)
    domain, geom = mesh.rectilinear(grid)
    domain = domain.withboundary(inflow="left", outflow="right", wall="top,bottom")

    # UNDER CONSTRUCTION

    # cloud of Gauss points
    gauss = domain.sample("gauss", degree=2)

    # Nutils namespace
    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = domain.basis("std", degree=1)  # linear finite elements
    ns.u = "basis_n ?lhs_n"  # solution
    ns.dudt = "basis_n (?lhs_n - ?lhs0_n) / ?dt"  # time derivative
    ns.vbasis = gauss.basis()
    ns.velocity_i = "vbasis_n ?velocity_ni"
    ns.k = 0.1  # diffusivity
    ns.xblob = 1, 1
    ns.uinit = ".5 - .5 tanh(((x_i - xblob_i) (x_i - xblob_i) - .5) / .1)"  # blob

    # define the weak form
    res = gauss.integral("(basis_n (dudt + (velocity_i u)_,i) + k basis_n,i u_,i) d:x" @ ns)

    # define Dirichlet boundary condition
    sqr = domain.boundary["inflow"].integral("u^2 d:x" @ ns, degree=2)
    cons = solver.optimize("lhs", sqr, droptol=1e-15)

    # preCICE setup
    interface = precice.Interface("Transport", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Transport-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    vertices = gauss.eval(ns.x)
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

    # coupling data
    velocity_id = interface.get_data_id("Velocity", mesh_id)

    precice_dt = interface.initialize()

    timestep = 0
    dt = 0.005

    # set blob as initial condition
    sqr = domain.integral("(u - uinit)^2" @ ns, degree=2)
    lhs0 = solver.optimize("lhs", sqr)

    # initialize the velocity values
    velocity_values = np.zeros_like(vertices)

    while interface.is_coupling_ongoing():

        if timestep % 1 == 0:  # visualize
            bezier = domain.sample("bezier", 2)
            x, u = bezier.eval(["x_i", "u"] @ ns, lhs=lhs0)
            with log.add(log.DataLog()):
                export.vtk("Transport_" + str(timestep), bezier.tri, x, T=u)

        # read velocity values from interface
        if interface.is_read_data_available():
            velocity_values = interface.read_block_vector_data(velocity_id, vertex_ids)

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)

        # solve nutils timestep
        lhs = solver.solve_linear(
            "lhs", res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt, velocity=velocity_values)
        )

        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        lhs0 = lhs

    interface.finalize()


if __name__ == "__main__":
    cli.run(main)
