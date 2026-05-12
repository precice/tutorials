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

    # COMMENT: nothing to change here
    log.info("Running Nutils")

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

    gauss = domain.sample("gauss", degree=4)

    ns = function.Namespace()
    ns.x = geom

    # COMMENT: looks like the weak is being defined here: the weak form now needs to incorporate the volume fraction
    # Note that preCICE/MercuryDPM provides a field called "Alpha", which is the solid volume fraction. The formula
    # 7.1 and 7.2 from the PDF use the fluid volume frection (epsilon_f), which is epsilon_f = (1 - Alpha).
    # Note that reading from preCICE makes only sense after the initialize call.
    # Moreover, we need to incorporate gravity (ε_f*ρ_f*g) in 7.2 and the drag force. The drag force is here zero either
    # Note that the drag force is called "ExplicitMomentum" in preCICE. I
    # would maybe later still try to implement the explicit-implicit split,
    # that's the only reason
    ns.ubasis = domain.basis("std", degree=2).vector(2)
    ns.pbasis = domain.basis("std", degree=1)
    ns.εbasis = ns.pbasis
    ns.αbasis = gauss.basis()
    ns.Fbasis = ns.αbasis.vector(2)
    ns.u_i = "ubasis_ni ?u_n"  # velocity
    ns.v_i = "ubasis_ni ?v_n"  # test velocity
    ns.u0_i = "ubasis_ni ?u0_n"  # previous velocity
    ns.p = "pbasis_n ?p_n"  # pressure
    ns.q = "pbasis_n ?q_n"  # test pressure
    ns.α = "αbasis_n ?α_n"  # solid fraction
    ns.ε = "εbasis_n ?ε_n"  # void fraction
    ns.ε0 = "εbasis_n ?ε0_n"  # previous void fraction
    ns.F_i = "Fbasis_ni ?F_n"  # drag force
    ns.g = np.array([0, -9.81])
    ns.ρf = 1.2  # fluid density
    ns.DεDt = "(ε - ε0) / ?dt + (ε u_i)_,i"
    ns.DuDt_i = "(u_i - u0_i) / ?dt + u_i,j u_j"
    ns.μ = 1.8e-5  # dynamic viscosity
    ns.uin = "10 x_1 (2 - x_1)"  # inflow profile

    # define the weak form
    # Gertjan implemented DεuDt_i instead of the ε DuDt_i
    # where ns.DεuDt_i = "(ε u_i - ε0 u0_i) / ?dt + ε u_i,j u_j"
    res = gauss.integral("(ε DuDt_i v_i - p (ε v_i)_,i / ρf + F_i v_i / ρf) d:x" @ ns)

    # the version closest to the paper notation
    # res += gauss.integral("((ε μ u_i,kk + μ u_i,j ε_,j) v_i) / ρf d:x" @ ns)
    # or Gertjan's equivalent version
    # res += gauss.integral("(μ (ε u_i,j)_,j v_i) / ρf d:x" @ ns)
    # But mistake in the paper? I think there shouldn't be the second derivative on u
    # Instead the standard integrated-by-parts form and a term with grad epsilon, i.e., from the product rule
    res += gauss.integral("(ε μ u_i,k v_i,k / ρf + (μ u_i,j ε_,j) v_i / ρf) d:x" @ ns)

    # Define the stabilization parameter \tau u according to Tezduyar, Eq. (16)
    ns.h = 0.08     # TODO: How do I get the element size here? note that this affects gamma
    ns.hconv = "h"  # element size related to convection (I would stick to h unless there is a better idea)
    ns.hdiff = "h"  # element size related to diffusion (I would stick to h here)
    ns.τu = "((1 / ?dt)^2 + (2 (u_k u_k + 1e-15)^0.5 / hconv)^2 + 9 (4 (μ / ρf) / hdiff^2)^2)^-0.5"

    # TODO: Refactor the \tau_u terms to use the same strong_res (remove duplication)
    # ns.strong_res_i = "ε DuDt_i + (ε p_,i / ρf) - (ε μ u_i,kk) + (F_i / ρf)"
    res += gauss.integral("(ε DuDt_i + ε p_,i / ρf - ε μ u_i,kk / ρf + F_i / ρf) τu v_i,j u_j d:x" @ ns)

    ns.γ = "μ / ρf + h (u_k u_k + 1e-15)^0.5"
    res += gauss.integral("γ DεDt v_i,i d:x" @ ns)

    # Continuity equation including the PSPG stabilization
    res += gauss.integral("q DεDt d:x" @ ns)
    res += gauss.integral("(ε DuDt_i + ε p_,i / ρf - ε μ u_i,kk / ρf + F_i / ρf) τu q_,i d:x" @ ns)

    # Dirichlet boundary condition
    sqr = domain.boundary["inflow"].integral("(u_0 - uin)^2 d:x" @ ns, degree=2)
    sqr += domain.boundary["inflow,outflow"].integral("u_1^2 d:x" @ ns, degree=2)
    sqr += domain.boundary["wall"].integral("u_k u_k d:x" @ ns, degree=2)
    cons = solver.optimize(["u"], sqr, droptol=1e-15)

    # preCICE setup
    participant = precice.Participant("Fluid", "../precice-config.xml", 0, 1)

    # define coupling mesh (nothing to change)
    mesh_name = "Fluid-Mesh"
    vertices = gauss.eval(ns.x)
    vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

    # COMMENT: names of the new coupling data fields in the preCICE configuration. The velocity was already there before
    solid_fraction_name = "Alpha"
    drag_force_name = "ExplicitMomentum"
    data_name = "Velocity"
    pressure_gradient_name = "PressureGradientFull"

    participant.initialize()

    timestep = 0
    solver_dt = 0.005
    precice_dt = participant.get_max_time_step_size()
    dt = min(precice_dt, solver_dt)

    # TODO: INITIAL CONDITION
    state = {
        'u': np.zeros(len(ns.ubasis)),
        'p': np.zeros(len(ns.pbasis)),  # for plotting
    }

    # Projection of the volume fraction
    sqr = gauss.integral('(ε - (1 - α))^2 d:x' @ ns)
    sys_project_ε = solver.System(sqr, trial='ε')

    sys_vans = solver.System(res, trial='u,p', test='v,q')

    bezier = domain.sample("bezier", 2)

    while participant.is_coupling_ongoing():

        if timestep > 1:  # visualize, in earlier timesteps, epsilon is not yet available
            x, u, p, eps = bezier.eval(["x_i", "u_i", "p", "ε"] @ ns, arguments=state)
            with log.add(log.DataLog()):
                export.vtk("Fluid_" + str(timestep), bezier.tri, x, u=u, p=p, eps=eps)

        precice_dt = participant.get_max_time_step_size()

        # potentially adjust non-matching timestep sizes
        dt = min(solver_dt, precice_dt)

        # COMMENT: here we would need to read the solid volume fraction and the drag force from preCICE and incorporate them into the weak form. The pressure gradient would also need to be written to preCICE, but let's start with the velocity first
        # Note that the read_data call takes a dt argument, which is the relative read time. For Euler backwards, it would simply be dt
        # solid_fraction_values = participant.read_data(mesh_name, solid_fraction_name, vertex_ids, ...)
        # drag_force_name = participant.read_data(mesh_name, drag_force_name, vertex_ids, ...)
        # Checkpointing for implicit coupling is generally not required

        state['α'] = participant.read_data(mesh_name, solid_fraction_name, vertex_ids, dt)
        state['F'] = participant.read_data(mesh_name, drag_force_name, vertex_ids, dt).flatten()

        # determine void ratio
        state = sys_project_ε.solve(arguments=state)

        # solve timestep
        state["u0"] = state["u"]
        state["ε0"] = state["ε"]
        state["dt"] = dt
        state = sys_vans.solve(constrain=cons, arguments=state, tol=1e-10)

        # COMMENT: here we would also need to write the pressure gradient to preCICE
        # Important note: OpenFOAM scales the pressue (and thereby also its gradient) by the fluid density, so we would
        # need to do the same scaling here. I set the fluid density to 1.2
        velocity_values, pressure_gradient_values = gauss.eval(['u_i', 'p_,i / ρf'] @ ns, arguments=state)
        participant.write_data(mesh_name, data_name, vertex_ids, velocity_values)
        participant.write_data(mesh_name, pressure_gradient_name, vertex_ids, pressure_gradient_values)

        # do the coupling
        participant.advance(dt)

        # advance variables
        timestep += 1

    participant.finalize()


if __name__ == "__main__":
    cli.run(main)
