# Turek benchmark fluid solver, modified from https://examples.nutils.org/official-turek/
#
# Thanks to the Technical University of Eindhoven for funding the development
# of this script.

from precice import Participant
from nutils import cli, export, function
from nutils.mesh import gmsh
from nutils.solver import System
from nutils.SI import Length, Density, Viscosity, Velocity, Time, Stress, Acceleration
from nutils.expression_v2 import Namespace
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
from itertools import count
import treelog as log
import numpy


# reference lenghts used in communication
precice_length = Length('m')
precice_time = Time('s')
precice_stress = Stress('Pa')


@dataclass
class Domain:
    channel_length: Length = Length('2.5m')
    channel_height: Length = Length('.41m')
    x_center: Length = Length('.2m')
    y_center: Length = Length('.2m')
    cylinder_radius: Length = Length('5cm')
    structure_length: Length = Length('35cm')
    structure_thickness: Length = Length('2cm')
    min_elemsize: Length = Length('4mm')
    max_elemsize: Length = Length('4cm')

    def generate_mesh(self):
        u = Length('m')  # temporary reference length

        topo, geom = gmsh(Path(__file__).parent / 'fluid.geo', dimension=2, order=2, numbers={
            'channel_length': self.channel_length / u,
            'channel_height': self.channel_height / u,
            'x_center': self.x_center / u,
            'y_center': self.y_center / u,
            'cylinder_radius': self.cylinder_radius / u,
            'structure_length': self.structure_length / u,
            'structure_thickness': self.structure_thickness / u,
            'min_elemsize': self.min_elemsize / u,
            'max_elemsize': self.max_elemsize / u})

        # consistency check
        # zeros = topo.boundary['inlet,wall,outlet,cylinder,structure'].integrate(function.normal(geom) * function.J(geom, 1), degree=2)
        # numpy.testing.assert_allclose(zeros, 0, atol=1e-14, err_msg='boundaries do not form a watertight hull')

        bezier = topo.sample('bezier', 2)
        bezier_structure = topo['fluid'].boundary['structure'].sample(
            'bezier', 3)
        bezier_cylinder = topo['fluid'].boundary['cylinder'].sample(
            'bezier', 3)
        A = topo.points['A'].sample('gauss', 1).eval(geom)
        with export.mplfigure('mesh.jpg', dpi=150) as fig:
            ax = fig.add_subplot(111)
            export.triplot(ax, bezier.eval(geom), hull=bezier.hull)
            export.triplot(ax, bezier_structure.eval(
                geom), hull=bezier_structure.tri, linewidth=1, linecolor='r')
            export.triplot(ax, bezier_cylinder.eval(
                geom), hull=bezier_cylinder.tri, linewidth=1, linecolor='b')
            ax.set_xlim(0, 2 * self.channel_height / u)

        return topo, geom * u


@dataclass
class Fluid:
    density: Density = Density('1kg/L')
    viscosity: Viscosity = Viscosity('1Pa*s')
    velocity: Velocity = Velocity('2m/s')

    def reynolds(self, reference_length):
        return self.density * self.velocity * reference_length / self.viscosity


@dataclass
class Dynamic:
    init: Time = Time('2s')
    window: Time = Time('1s')
    gamma: float = .8
    beta: float = .4

    def __post_init__(self):
        self.timeseries = defaultdict(deque)

    def ramp_up(self, t):
        return .5 - .5 * numpy.cos(numpy.pi * min(t / self.init, 1))

    def add_and_plot(self, name, t, tunit, v, vunit, ax):
        'Add data point and plot time series for past window.'

        d = self.timeseries[name]
        d.append((t, v))
        while t - d[0][0] > self.window:
            d.popleft()
        times = numpy.array([t / tunit for t, _ in d])
        values = numpy.array([v / vunit for _, v in d])
        ax.plot(times, values)
        ax.set_ylabel(f'{name} [{vunit}]')
        ax.grid()
        ax.autoscale(enable=True, axis='x', tight=True)

    # The Newmark-beta scheme is used for time integration, taking displacement
    # (solid) or velocity (fluid) as primary variable with time derivatives
    # introduced via helper arguments that are updated after every solve.
    #
    # d = d0 + δt u0 + .5 δt^2 aβ, where aβ = (1-2β) a0 + 2β a
    # => δd = δt u0 + δt^2 [ .5 a0 + β δa ]
    # => δa = [ δd / δt^2 - u0 / δt - .5 a0 ] / β
    #
    # u = u0 + δt aγ, where aγ = (1-γ) a0 + γ a
    # => δu = δt [ a0 + γ δa ]
    # => δa = [ δu / δt - a0 ] / γ

    def newmark_defo(self, δt, δd, u0, a0):
        'Compute first and second time derivative of deformation.'
        δa = (δd / δt**2 - u0 / δt - .5 * a0) / self.beta
        δu = δt * (a0 + self.gamma * δa)
        return u0 + δu, a0 + δa

    def newmark_defo_fields(self, d):
        'Generate velocity and acceleration field.'
        # By using the same scaling factor for time and time derivatives we can
        # apply newmark_defo directly on the non-dimensional coefficients in
        # newmark_update.
        return self.newmark_defo(
            (function.field('t') - function.field('t0')) * precice_time,
            d - function.replace_arguments(d, 'd:d0'),
            function.replace_arguments(d, 'd:ḋ0') / precice_time,
            function.replace_arguments(d, 'd:d̈0') / precice_time**2)

    def newmark_velo(self, δt, δu, a0):
        'Compute first time derivative of velocity.'
        δa = (δu / δt - a0) / self.gamma
        return a0 + δa

    def newmark_velo_field(self, u):
        'Generate acceleration field.'
        # By using the same scaling factor for time and time derivatives we can
        # apply newmark_velo directly on the non-dimensional coefficients in
        # newmark_update.
        return self.newmark_velo(
            (function.field('t') - function.field('t0')) * precice_time,
            u - function.replace_arguments(u, 'u:u0'),
            function.replace_arguments(u, 'u:u̇0') / precice_time)

    def newmark_update(self, δt, t, d, u, t0=-1., d0=0., ḋ0=0., d̈0=0., u0=0., u̇0=0., **args):
        'Construct time derivatives, update time, construct initial guess.'
        # Note: all arguments are non-dimensional.
        ḋ, d̈ = self.newmark_defo(t - t0, d - d0, ḋ0, d̈0)
        u̇ = self.newmark_velo(t - t0, u - u0, u̇0)
        d_predict = d + δt * ḋ + .5 * δt**2 * d̈
        u_predict = u + δt * u̇
        return dict(args, t0=t, t=t + δt, d=d_predict, d0=d, ḋ0=ḋ, d̈0=d̈, u=u_predict, u0=u, u̇0=u̇)


def main(domain: Domain = Domain(), fluid: Fluid = Fluid(), dynamic: Dynamic = Dynamic()):

    participant = Participant('Fluid', '../precice-config.xml', 0, 1)

    log.info('Re:', fluid.reynolds(2 * domain.cylinder_radius))

    topo, geom = domain.generate_mesh()

    d = topo.field('d', btype='std', degree=1, shape=[2]) * precice_length

    sqr = topo.boundary['structure'].sample(
        'bezier', 2).integral((d - geom) @ (d - geom)) / 'm2'
    r_points = System(sqr, trial='d').solve_constraints(droptol=1e-10)['d']
    r_where = numpy.isfinite(r_points).any(1)
    assert not numpy.isnan(r_points[r_where]).any()
    r_name = "Fluid-Mesh-Nodes"
    r_ids = participant.set_mesh_vertices(r_name, r_points[r_where])

    w_name = "Fluid-Mesh-Centers"
    w_sample = topo.boundary['structure'].sample('gauss', degree=1)
    w_ids = participant.set_mesh_vertices(
        w_name, w_sample.eval(geom) / precice_length)

    participant.initialize()

    ns = Namespace()
    ns.δ = function.eye(2)
    ns.ρf = fluid.density
    ns.μf = fluid.viscosity

    ns.xref = geom
    ns.define_for('xref', gradient='∇ref', jacobians=('dVref', 'dSref'))

    ns.d = d
    ns.x_i = 'xref_i + d_i'
    ns.F_ij = '∇ref_j(x_i)'  # deformation gradient tensor
    ns.C_ij = 'F_ki F_kj'  # right Cauchy-Green deformation tensor
    ns.J = numpy.linalg.det(ns.F)

    ns.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))

    ns.urel = topo.field('u', btype='std', degree=2,
                         shape=(2,)) * fluid.velocity
    ns.v, ns.a = dynamic.newmark_defo_fields(ns.d)
    ns.arel = dynamic.newmark_velo_field(ns.urel)
    ns.u_i = 'v_i + urel_i'
    ns.DuDt_i = 'a_i + arel_i + ∇_j(u_i) urel_j'  # material derivative

    ns.p = topo.field('p', btype='std', degree=1) * \
        fluid.viscosity * fluid.velocity / domain.cylinder_radius
    ns.σ_ij = 'μf (∇_j(u_i) + ∇_i(u_j)) - p δ_ij'  # fluid stress tensor

    y = ns.xref[1] / domain.channel_height
    uin = 6 * fluid.velocity * y * (1 - y)
    sqr = topo.boundary['wall,cylinder,structure'].integral(
        ns.urel @ ns.urel, degree=4) / 'm2/s2'
    sqr += topo.boundary['wall,cylinder,inlet,outlet'].integral(
        ns.d @ ns.d, degree=4) / 'm2'
    sqr += topo.boundary['inlet'].integral(
        (ns.urel[0] - uin)**2 * ns.dSref, degree=4) / 'm3/s2'
    sqr += topo.boundary['inlet,outlet'].integral(
        ns.urel[1]**2, degree=4) / 'm2/s2'
    cons = System(sqr, trial='u,d').solve_constraints(droptol=1e-10)
    ucons = cons['u']
    cons['u'] = ucons * 0

    # ρf DuDt = div σ =>  ∀q: 0 = ∫ q·(ρf DuDt - div σ) = ∫ (ρf q·DuDt + ∇q:σ)
    ns.utest = function.replace_arguments(
        ns.urel, 'u:utest') / (fluid.viscosity * fluid.velocity**2)
    ns.ptest = function.replace_arguments(
        ns.p, 'p:ptest') / (fluid.viscosity * fluid.velocity**2)
    res = topo.integral(
        '(utest_i ρf DuDt_i + ∇_j(utest_i) σ_ij + ptest ∇_k(u_k)) dV' @ ns, degree=4)
    sys_up = System(res, trial='u,p', test='utest,ptest')

    # t = -σ·n => ∀q: ∮ -q·t = ∮ q·σ·n = ∫ div(q·σ) = ∫ (∇q:σ + q·div σ) = ∫ (∇q:σ + ρf q·DuDt)
    ns.traction = topo.field('traction', btype='std', degree=2, shape=[
        2]) * (fluid.viscosity * fluid.velocity / domain.cylinder_radius)
    res += topo.boundary['cylinder,structure'].integral(
        'utest_i traction_i dS' @ ns, degree=4)
    sys_t = System(res, trial='traction', test='utest')
    sqr = topo.boundary['cylinder,structure'].integral(
        'traction_i traction_i' @ ns, degree=4) / 'Pa2'
    cons['traction'] = numpy.isnan(
        System(sqr, trial='traction').solve_constraints(droptol=1e-10)['traction'])
    F = topo.boundary['cylinder,structure'].integral('traction_i dS' @ ns, degree=4)

    # mesh continuation with jacobian based stiffness
    sqr = topo.integral('C_kk - 2 log(J)' @ ns, degree=4)
    sys_d = System(sqr, trial='d')

    # initial values
    args = {a: numpy.zeros(function.arguments_for(res)[a].shape) for a in 'udt'}

    bezier = topo['fluid'].sample('bezier', 3)
    bezier = bezier.subset(bezier.eval(geom[0]) < 2.2 * domain.channel_height)
    x_bz = bezier.bind(ns.x)
    u_bz = bezier.bind(ns.u)
    p_bz = bezier.bind(ns.p) - \
        topo.points['B'].sample('gauss', 1).bind(ns.p)[0]

    bbezier = topo['fluid'].boundary['cylinder,structure'].sample('bezier', 3)
    x_bbz = bbezier.bind(ns.x)

    with log.iter.plain('time step', count()) as iter_context:
        while participant.is_coupling_ongoing():

            if participant.requires_writing_checkpoint():
                checkpoint = args

            timestep = participant.get_max_time_step_size()

            # receive flap displacements
            cons['d'][r_where] = participant.read_data(
                r_name, 'Displacement', r_ids, timestep)
            # generate inflow constraints
            cons['u'] = ucons * dynamic.ramp_up(args['t'] * precice_time)
            # update time derivatives (implies copy)
            args = dynamic.newmark_update(timestep, **args)
            # extend mesh
            args = sys_d.solve(arguments=args, constrain=cons, tol=1e-10)
            # solve Navier-Stokes time step
            args = sys_up.solve(arguments=args, constrain=cons, tol=1e-10)
            # project boundary tractions
            args = sys_t.solve(arguments=args, constrain=cons, tol=1e-10)
            # send boundary tractions
            participant.write_data(
                w_name, 'Stress', w_ids, w_sample.eval(ns.traction, args) / precice_stress)

            participant.advance(timestep)

            if participant.requires_reading_checkpoint():
                args = checkpoint

            if participant.is_time_window_complete():
                next(iter_context)

                t = args['t'] * precice_time

                x, xb, u, p = function.eval([x_bz, x_bbz, u_bz, p_bz], args)
                with export.mplfigure('solution.jpg', dpi=150) as fig:
                    pstep = 25 * fluid.viscosity * fluid.velocity / domain.channel_height
                    ax = fig.add_subplot(
                        111, title=f'flow at t={t:.3s}, pressure contours every {pstep:.0Pa}', ylabel='[m]')
                    vmax = 2 * fluid.velocity * dynamic.ramp_up(t)
                    im = export.triplot(ax, x / 'm', numpy.linalg.norm(u / 'm/s', axis=1),
                                        tri=bezier.tri, cmap='inferno', clim=(0, vmax / 'm/s'))
                    ax.tricontour(*(x / 'm').T, bezier.tri, p / pstep, numpy.arange(*numpy.quantile(numpy.ceil(p / pstep), [0, 1])),
                                  colors='white', linestyles='solid', linewidths=1, alpha=.33)
                    fig.colorbar(im, orientation='horizontal',
                                 label=f'velocity [m/s]')
                    export.triplot(ax, xb / 'm', hull=bbezier.tri, linewidth=1)
                    ax.set_xlim(0, 2 * domain.channel_height / 'm')
                    ax.set_ylim(0, domain.channel_height / 'm')

                D, L = function.eval(F, arguments=args)
                log.info(f'lift: {L:N/m}')
                log.info(f'drag: {D:N/m}')
                with export.mplfigure('force.jpg', dpi=150) as fig:
                    dynamic.add_and_plot(
                        'lift', t, 's', L, 'N/m', ax=fig.add_subplot(211))
                    dynamic.add_and_plot(
                        'drag', t, 's', D, 'N/m', ax=fig.add_subplot(212, xlabel='time [s]'))

    participant.finalize()


if __name__ == '__main__':
    cli.run(main)
