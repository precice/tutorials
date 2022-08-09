#! /usr/bin/env python3

#
# Advection-Diffusion equation for a single species with a velocity field read from preCICE on the complete volume.
#

from nutils import function, mesh, cli, solver, export
import treelog as log
import numpy as np
import precice


def reinitialize_namespace(domain, geom):
    # cloud of Gauss points
    gauss = domain.sample("gauss", degree=2)

    # Nutils namespace
    ns = function.Namespace(fallback_length=2)
    ns.x = geom
    ns.basis = domain.basis("h-std", degree=1)  # linear finite elements
    ns.u = "basis_n ?solu_n"  # solution
    ns.projectedu = "basis_n ?projectedsolu_n"
    ns.gradu = "u_,i"  # gradient of solution
    ns.dudt = "basis_n (?solu_n - ?solu0_n) / ?dt"  # time derivative
    ns.vbasis = gauss.basis()
    ns.velocity_i = "vbasis_n ?velocity_ni"
    ns.k = 0.1  # diffusivity
    ns.xblob = 1, 1
    ns.uinit = ".5 - .5 tanh(((x_i - xblob_i) (x_i - xblob_i) - .5) / .1)"  # blob

    # define the weak form
    res = gauss.integral("(basis_n (dudt + (velocity_i u)_,i) + k basis_n,i u_,i) d:x" @ ns)

    # define Dirichlet boundary condition
    sqr = domain.boundary["inflow"].integral("u^2 d:x" @ ns, degree=2)
    cons = solver.optimize("solu", sqr, droptol=1e-15)

    return ns, res, cons, gauss


def refine_mesh(ns, domain_coarse, domain_nm1, solu_nm1):
    """
    At the time of the calling of this function a predicted solution exists in ns.phi
    """
    # ----- Refine the coarse mesh according to the projected solution to get a predicted refined topology ----
    domain_ref = domain_coarse
    for level in range(3):
        print("refinement level = {}".format(level))
        domain_union1 = domain_nm1 & domain_ref
        smpl = domain_union1.sample('uniform', 5)
        ielem, criterion = smpl.eval([domain_ref.f_index, function.sqrt(ns.gradu[0]**2 + ns.gradu[1]**2) > 2.0],
                                     solu=solu_nm1)

        # Refine the elements for which at least one point tests true.
        domain_ref = domain_ref.refined_by(np.unique(ielem[criterion]))
        # ----------------------------------------------------------------------------------------------------

    # Create a new projection mesh which is the union of the previous refined mesh and the predicted mesh
    domain_union = domain_nm1 & domain_ref

    # ----- Project the solution of the last time step on the projection mesh -----
    ns.projectedu = function.dotarg('projectedsolu', domain_ref.basis('h-std', degree=1))
    sqru = domain_union.integral((ns.projectedu - ns.u) ** 2, degree=2)
    solu = solver.optimize('projectedsolu', sqru, droptol=1E-12, arguments=dict(solu=solu_nm1))

    return domain_ref, solu


def main():

    print("Running Nutils")

    # Remeshing
    n_remeshing = 10

    # define the Nutils mesh
    nx = 60
    ny = 16
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
    ns.u = "basis_n ?solu_n"  # solution
    ns.projectedu = "basis_n ?projectedsolu_n"
    ns.gradu = "u_,i"  # gradient of solution
    ns.dudt = "basis_n (?solu_n - ?solu0_n) / ?dt"  # time derivative
    ns.vbasis = gauss.basis()
    ns.velocity_i = "vbasis_n ?velocity_ni"
    ns.k = 0.1  # diffusivity
    ns.xblob = 1, 1
    ns.uinit = ".5 - .5 tanh(((x_i - xblob_i) (x_i - xblob_i) - .5) / .1)"  # blob

    # define the weak form
    res = gauss.integral("(basis_n (dudt + (velocity_i u)_,i) + k basis_n,i u_,i) d:x" @ ns)

    # define Dirichlet boundary condition
    sqr = domain.boundary["inflow"].integral("u^2 d:x" @ ns, degree=2)
    cons = solver.optimize("solu", sqr, droptol=1e-15)

    timestep = 0
    dt = 0.005

    # set blob as initial condition
    sqr = domain.integral("(u - uinit)^2" @ ns, degree=2)
    solu0 = solver.optimize("solu", sqr)

    # Initial refinement according to initial condition
    print("Performing initial mesh refinement")
    for level in range(3):
        print("refinement level = {}".format(level))
        smpl = domain.sample('uniform', 5)
        ielem, criterion = smpl.eval([domain.f_index, function.sqrt(ns.gradu[0]**2 + ns.gradu[1]**2) > 2.0], solu=solu0)

        # Refine the elements for which at least one point tests true.
        domain = domain.refined_by(np.unique(ielem[criterion]))

        ns, res, cons, gauss = reinitialize_namespace(domain, geom)

        # set blob as initial condition after each refinement
        sqr = domain.integral("(u - uinit)^2" @ ns, degree=2)
        solu0 = solver.optimize("solu", sqr)

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

    # initialize the velocity values
    velocity_values = np.zeros_like(vertices)

    solu = solu0

    while interface.is_coupling_ongoing():
        if timestep % n_remeshing == 0:
            domain, solu = refine_mesh(ns, domain_coarse, domain, solu)
            ns, res, cons, gauss = reinitialize_namespace(domain, geom)

            vertices = gauss.eval(ns.x)
            interface.reset_mesh(mesh_id)  # Throws away the entire mesh
            vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)  # Redefine the mesh

        if timestep % 1 == 0:  # visualize
            bezier = domain.sample("bezier", 2)
            x, u = bezier.eval(["x_i", "u"] @ ns, solu=solu0)
            with log.add(log.DataLog()):
                export.vtk("Transport_" + str(timestep), bezier.tri, x, T=u)

        # read velocity values from interface
        if interface.is_read_data_available():
            velocity_values = interface.read_block_vector_data(velocity_id, vertex_ids)

        # potentially adjust non-matching timestep sizes
        dt = min(dt, precice_dt)

        # solve nutils timestep
        args = dict(solu0=solu0, dt=dt, velocity=velocity_values)
        solu = solver.solve_linear("solu", res, constrain=cons, arguments=args)

        # do the coupling
        precice_dt = interface.advance(dt)

        # advance variables
        timestep += 1
        solu0 = solu

    interface.finalize()


if __name__ == "__main__":
    cli.run(main)
