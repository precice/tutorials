#! /usr/bin/env python3

#
# Incompressible NSE solved within a channel geometry with parabolic inflow profile and an obstacle attached to the bottom towards the middle of the domain. The fluid field is initialized with a Stokes solution. The resulting velocity field is written to preCICE on the complete volume.
#

from nutils import function, mesh, cli, solver, export
import treelog as log
import numpy as np
import precice
from mpi4py import MPI


def main():

    print("Running Nutils")

    # define the Nutils mesh
    nx = 48
    ny = 16
    step_start = nx // 3
    step_end = nx // 2
    step_hight = ny // 2

    grid = np.linspace(0, 6, nx + 1), np.linspace(0, 2, ny + 1)
    domain, geom = mesh.rectilinear(grid)
    domain = domain.withboundary(inflow="left", outflow="right", wall="top,bottom") - domain[
        step_start:step_end, :step_hight
    ].withboundary(wall="left,top,right")

    # cloud of Gauss points
    gauss = domain.sample("gauss", degree=4)

    # Nutils namespace
    ns = function.Namespace()
    ns.x = geom

    ns.ubasis = domain.basis("std", degree=2).vector(2)
    ns.pbasis = domain.basis("std", degree=1)
    ns.u_i = "ubasis_ni ?u_n"  # solution
    ns.p = "pbasis_n ?p_n"  # solution
    ns.dudt_i = "ubasis_ni (?u_n - ?u0_n) / ?dt"  # time derivative
    ns.μ = 0.5  # viscosity
    ns.σ_ij = "μ (u_i,j + u_j,i) - p δ_ij"
    ns.uin = "10 x_1 (2 - x_1)"  # inflow profile

    # define the weak form, Stokes problem
    ures = gauss.integral("ubasis_ni,j σ_ij d:x" @ ns)
    pres = gauss.integral("pbasis_n u_k,k d:x" @ ns)

    # define Dirichlet boundary condition
    sqr = domain.boundary["inflow"].integral("(u_0 - uin)^2 d:x" @ ns, degree=2)
    sqr += domain.boundary["inflow,outflow"].integral("u_1^2 d:x" @ ns, degree=2)
    sqr += domain.boundary["wall"].integral("u_k u_k d:x" @ ns, degree=2)
    cons = solver.optimize(["u"], sqr, droptol=1e-15)

    # preCICE setup
    participant = precice.Participant("Fluid", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Fluid-Mesh"
    vertices = gauss.eval(ns.x)
    vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

    # coupling data
    data_name = "Velocity"

    participant.initialize()

    timestep = 0
    solver_dt = 0.005
    precice_dt = participant.get_max_time_step_size()
    dt = min(precice_dt, solver_dt)

    state = solver.solve_linear(("u", "p"), (ures, pres), constrain=cons)  # initial condition

    # add convective term and time derivative for Navier-Stokes
    ures += gauss.integral("ubasis_ni (dudt_i + μ (u_i u_j)_,j) d:x" @ ns)

    while participant.is_coupling_ongoing():

        if timestep % 1 == 0:  # visualize
            bezier = domain.sample("bezier", 2)
            x, u, p = bezier.eval(["x_i", "u_i", "p"] @ ns, **state)
            with log.add(log.DataLog()):
                export.vtk("Fluid_" + str(timestep), bezier.tri, x, u=u, p=p)

        precice_dt = participant.get_max_time_step_size()

        # potentially adjust non-matching timestep sizes
        dt = min(solver_dt, precice_dt)

        # solve Nutils timestep
        state["u0"] = state["u"]
        state["dt"] = dt
        state = solver.newton(("u", "p"), (ures, pres), constrain=cons, arguments=state).solve(1e-10)

        velocity_values = gauss.eval(ns.u, **state)
        participant.write_data(mesh_name, data_name, vertex_ids, velocity_values)

        # do the coupling
        participant.advance(dt)

        # advance variables
        timestep += 1

    participant.finalize()


if __name__ == "__main__":
    cli.run(main)
