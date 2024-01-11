#! /usr/bin/env python3

from nutils import mesh, function, solver, export, cli
from nutils.expression_v2 import Namespace
import numpy
import treelog
import precice


def main(young=4e6, density=3e3, poisson=0.3, nelems=2, solver_dt=0.01, npoints_per_elem=3):

    topo, geom = mesh.rectilinear([numpy.linspace(-0.05, 0.05, nelems + 1), numpy.linspace(0, 1, 10 * nelems + 1)])
    wall = topo.boundary['left,top,right'].sample('uniform', npoints_per_elem)

    ns = Namespace()
    ns.X = geom
    ns.δ = function.eye(2)
    ns.define_for('X', gradient='∇', normal='n', jacobians=('dV', 'dS'))
    ns.add_field('dt')
    ns.add_field(('u', 'u0', 'testu', 'v', 'v0', 'testv'), topo.basis('std', degree=1), shape=(2,))
    ns.add_field('F', wall.basis(), shape=(2,))
    ns.qw = 1 / npoints_per_elem  # quadrature weight
    ns.t_i = 'F_i / qw dS'
    ns.dudt_i = '(u_i - u0_i) / dt'
    ns.dvdt_i = '(v_i - v0_i) / dt'
    ns.ρ = density
    ns.λ = young * poisson / ((1 + poisson) * (1 - 2 * poisson))
    ns.μ = young / (2 * (1 + poisson))
    ns.σ_ij = 'λ ∇_k(u_k) δ_ij + μ (∇_j(u_i) + ∇_i(u_j))'

    # make sure we correctly scale point forces to tractions
    test_force = numpy.random.normal(size=(wall.npoints, 2))
    numpy.testing.assert_almost_equal(
        actual=wall.integrate('t_i dS' @ ns, F=test_force),
        desired=test_force.sum(0),
        decimal=10,
        err_msg='nutils error: failed to recover net force',
    )

    # continuum equations: ρ v' = ∇·σ + F, u' = v
    res = topo.integral('testv_i (dudt_i - v_i) dV' @ ns, degree=2)
    res += topo.integral('(testu_i ρ dvdt_i + ∇_j(testu_i) σ_ij) dV' @ ns, degree=2)
    res -= wall.integral('testu_i t_i dS' @ ns)

    # boundary conditions: fully constrained at y=0
    sqr = topo.boundary['bottom'].integral('u_k u_k' @ ns, degree=2)
    cons = solver.optimize('u,', sqr, droptol=1e-10)

    # initial conditions: undeformed and unmoving
    sqr = topo.integral('u_k u_k + v_k v_k' @ ns, degree=2)
    arguments = solver.optimize('u,v', sqr, constrain=cons)

    # preCICE setup
    solver_process_index = 0
    solver_process_size = 1
    participant = precice.Participant("Solid", "../precice-config.xml", solver_process_index, solver_process_size)
    mesh_name = "Solid-Mesh"
    vertex_ids = participant.set_mesh_vertices(mesh_name, wall.eval(ns.X))
    write_data_name = "Displacement"
    read_data_name = "Force"

    # initialize preCICE
    participant.initialize()

    timestep = 0
    force = numpy.zeros((wall.npoints, 2))

    while participant.is_coupling_ongoing():
        with treelog.context(f'timestep {timestep}'):

            # save checkpoint
            if participant.requires_writing_checkpoint():
                checkpoint = timestep, arguments

            precice_dt = participant.get_max_time_step_size()
            dt = min(precice_dt, solver_dt)

            # read forces from participant at the end of the timestep
            force = participant.read_data(mesh_name, read_data_name, vertex_ids, dt)

            # advance variables
            timestep += 1

            arguments = dict(dt=dt, u0=arguments['u'], v0=arguments['v'], F=force)
            arguments = solver.solve_linear('u:testu,v:testv', res, arguments=arguments, constrain=cons)

            # write displacements to participant
            write_data = wall.eval(ns.u, **arguments)
            participant.write_data(mesh_name, write_data_name, vertex_ids, write_data)

            # do the coupling
            participant.advance(dt)

            # read checkpoint if required
            if participant.requires_reading_checkpoint():
                timestep, arguments = checkpoint

    participant.finalize()


if __name__ == '__main__':
    cli.run(main)
