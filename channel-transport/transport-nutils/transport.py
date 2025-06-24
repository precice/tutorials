#! /usr/bin/env python3

#
# Advection-Diffusion equation for a single species with a velocity field read from preCICE on the complete volume.
#

from nutils import function, mesh, cli, solver, export
import treelog as log
import argparse
import numpy as np
import json
import precice
from mpi4py import MPI


def reinitialize_namespace(domain, geom):
    # cloud of Gauss points
    gauss = domain.sample("bezier", degree=2)

    # Nutils namespace
    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = domain.basis("h-std", degree=1)  # linear finite elements
    ns.u = "basis_n ?lhs_n"  # solution
    ns.projectedu = "basis_n ?projectedlhs_n"
    ns.gradu = "u_,i"  # gradient of lhstion
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

    return ns, res, cons, gauss


def refine_mesh(limits, ns, domain_coarse, domain_nm1, lhs_nm1):
    """
    At the time of the calling of this function a predicted lhstion exists in ns.phi
    """
    # ----- Refine the coarse mesh according to the projected lhstion to get a predicted refined topology ----
    domain_ref = domain_coarse
    for level, limit in enumerate(limits):
        print("refinement level = {}".format(level))
        domain_union1 = domain_nm1 & domain_ref
        smpl = domain_union1.sample('uniform', 5)
        ielem, criterion = smpl.eval([domain_ref.f_index, function.sqrt(ns.gradu[0]**2 + ns.gradu[1]**2) > limit],
                                     lhs=lhs_nm1)

        # Refine the elements for which at least one point tests true.
        domain_ref = domain_ref.refined_by(np.unique(ielem[criterion]))
        # ----------------------------------------------------------------------------------------------------

    # Create a new projection mesh which is the union of the previous refined mesh and the predicted mesh
    domain_union = domain_nm1 & domain_ref

    # ----- Project the lhstion of the last time step on the projection mesh -----
    ns.projectedu = function.dotarg('projectedlhs', domain_ref.basis('h-std', degree=1))
    sqru = domain_union.integral((ns.projectedu - ns.u) ** 2, degree=2)
    lhs = solver.optimize('projectedlhs', sqru, droptol=1E-12, arguments=dict(lhs=lhs_nm1))

    return domain_ref, lhs


def main(remesh: bool, frequency: int, limits: [float], visualize: bool):

    print("Running Nutils")

    if frequency > 0:
        # Use a coarser grid for remeshing cases
        nx = 60
        ny = 16
    else:
        # define the Nutils mesh
        nx = 120
        ny = 32

    step_start = nx // 3
    step_end = nx // 2
    step_height = ny // 2

    grid = np.linspace(0, 6, nx + 1), np.linspace(0, 2, ny + 1)
    domain, geom = mesh.rectilinear(grid)
    domain = domain.withboundary(inflow="left", outflow="right", wall="top,bottom") - domain[
        step_start:step_end, :step_height
    ].withboundary(wall="left,top,right")
    domain_coarse = domain  # Retain the original coarse domain for mesh refinement later on

    # cloud of Gauss points
    gauss = domain.sample("gauss", degree=2)

    # Nutils namespace
    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = domain.basis("std", degree=1)  # linear finite elements
    ns.u = "basis_n ?lhs_n"  # solution
    ns.projectedu = "basis_n ?projectedlhs_n"
    ns.gradu = "u_,i"  # gradient of lhstion
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

    # set blob as initial condition
    sqr = domain.integral("(u - uinit)^2" @ ns, degree=2)
    lhs0 = solver.optimize("lhs", sqr)

    if remesh:
        # Initial refinement according to initial condition
        print("Performing initial mesh refinement")
        for level, limit in enumerate(limits):
            print("refinement level = {}".format(level))
            smpl = domain.sample('uniform', 5)
            ielem, criterion = smpl.eval([domain.f_index, function.sqrt(
                ns.gradu[0]**2 + ns.gradu[1]**2) > limit], lhs=lhs0)

            # Refine the elements for which at least one point tests true.
            domain = domain.refined_by(np.unique(ielem[criterion]))

            ns, res, cons, gauss = reinitialize_namespace(domain, geom)

            # set blob as initial condition after each refinement
            sqr = domain.integral("(u - uinit)^2" @ ns, degree=2)
            lhs0 = solver.optimize("lhs", sqr)

    lhs = lhs0

    # preCICE setup
    participant = precice.Participant("Transport", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Transport-Mesh"
    vertices = gauss.eval(ns.x)
    vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

    # coupling data
    data_name = "Velocity"

    participant.initialize()

    timestep = 1
    solver_dt = 0.005
    precice_dt = participant.get_max_time_step_size()
    dt = min(precice_dt, solver_dt)
    t = 0

    # initialize the velocity values
    velocity_values = np.zeros_like(vertices)

    series = {"file-series-version": "1.0", "files": []}

    while participant.is_coupling_ongoing():

        precice_dt = participant.get_max_time_step_size()

        # potentially adjust non-matching timestep sizes
        dt = min(solver_dt, precice_dt)

        # read velocity values from participant
        velocity_values = participant.read_data(mesh_name, data_name, vertex_ids, dt)

        # solve nutils timestep
        lhs = solver.solve_linear(
            "lhs", res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt, velocity=velocity_values)
        )

        if remesh and timestep % frequency == 0:
            domain, lhs = refine_mesh(limits, ns, domain_coarse, domain, lhs)
            ns, res, cons, gauss = reinitialize_namespace(domain, geom)

            vertices = gauss.eval(ns.x)
            participant.reset_mesh(mesh_name)  # Throws away the entire mesh
            vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)  # Redefine the mesh

        if visualize:
            bezier = domain.sample("bezier", 2)
            x, u = bezier.eval(["x_i", "u"] @ ns, lhs=lhs)
            with log.add(log.DataLog()):
                export.vtk(f"Transport_{timestep}", bezier.tri, x, T=u)
                series["files"].append({"name": f"Transport_{timestep}.vtk", "time": t})

        # do the coupling
        participant.advance(dt)

        # advance variables
        t += dt
        timestep += 1
        lhs0 = lhs

    if visualize:
        with open("Transport.vtk.series", "w") as f:
            json.dump(series, f)

    participant.finalize()


if __name__ == "__main__":

    def run(remesh: bool = False, frequency: int = 2, visualize: bool = True):
        # we use 1 and 2 as default refinement limits.
        main(remesh, frequency, [1.0, 2.0], visualize)

    cli.run(run)
