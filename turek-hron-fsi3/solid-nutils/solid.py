# Turek benchmark solid solver, modified from https://examples.nutils.org/official-turek/
#
# Thanks to the Technical University of Eindhoven for funding the development
# of this script.

from precice import Participant
from nutils import cli, export, function
from nutils.mesh import gmsh
from nutils.solver import System
from nutils.SI import Length, Density, Viscosity, Velocity, Time, Stress, Acceleration, Pressure
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
    x_center: Length = Length('.2m')
    y_center: Length = Length('.2m')
    cylinder_radius: Length = Length('5cm')
    structure_length: Length = Length('35cm')
    structure_thickness: Length = Length('2cm')
    elemsize: Length = Length('4mm')

    def generate_mesh(self):
        u = Length('m')  # temporary reference length

        topo, geom = gmsh(Path(__file__).parent / 'solid.geo', dimension=2, order=2, numbers={
            'x_center': self.x_center / u,
            'y_center': self.y_center / u,
            'cylinder_radius': self.cylinder_radius / u,
            'structure_length': self.structure_length / u,
            'structure_thickness': self.structure_thickness / u,
            'elemsize': self.elemsize / u})

        bezier = topo.sample('bezier', 2)
        bezier_fluid = topo.boundary['fluid'].sample('bezier', 3)
        bezier_cylinder = topo.boundary['cylinder'].sample('bezier', 3)
        with export.mplfigure('mesh.jpg', dpi=150) as fig:
            ax = fig.add_subplot(111)
            export.triplot(ax, bezier.eval(geom), hull=bezier.hull)
            export.triplot(ax, bezier_fluid.eval(geom),
                           hull=bezier_fluid.tri, linewidth=1, linecolor='r')
            export.triplot(ax, bezier_cylinder.eval(
                geom), hull=bezier_cylinder.tri, linewidth=1, linecolor='b')
            ax.set_xlim(0, 4 * self.x_center / u)
            ax.set_ylim(0, 2 * self.y_center / u)

        return topo, geom * u


@dataclass
class Solid:
    density: Density = Density('1kg/L')
    poisson_ratio: float = .4
    shear_modulus: Pressure = Pressure('2MPa')
    gravity: Acceleration = Acceleration('0m/s2')

    def lame_parameters(self):
        return 2 * self.shear_modulus * self.poisson_ratio / (1 - 2 * self.poisson_ratio), self.shear_modulus

    def young(self):
        return 2 * self.shear_modulus * (1 + self.poisson_ratio)


@dataclass
class Dynamic:
    window: Time = Time('1s')
    gamma: float = .8
    beta: float = .4

    def __post_init__(self):
        self.timeseries = defaultdict(deque)

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

    def newmark_update(self, δt, t, d, t0=-1., d0=0., ḋ0=0., d̈0=0., **args):
        'Construct time derivatives, update time, construct initial guess.'
        # Note: all arguments are non-dimensional.
        ḋ, d̈ = self.newmark_defo(t - t0, d - d0, ḋ0, d̈0)
        d_predict = d + δt * ḋ + .5 * δt**2 * d̈
        return dict(args, t0=t, t=t + δt, d=d_predict, d0=d, ḋ0=ḋ, d̈0=d̈)


def main(domain: Domain = Domain(), solid: Solid = Solid(), dynamic: Dynamic = Dynamic()):

    participant = Participant('Solid', '../precice-config.xml', 0, 1)

    topo, geom = domain.generate_mesh()

    rw_name = "Solid-Mesh"
    rw_sample = topo.boundary.sample('gauss', degree=1)
    rw_ids = participant.set_mesh_vertices(
        rw_name, rw_sample.eval(geom) / precice_length)

    participant.initialize()

    ns = Namespace()
    ns.δ = function.eye(2)
    ns.xref = geom
    ns.define_for('xref', gradient='∇ref', jacobians=('dVref', 'dSref'))

    ns.ρs = solid.density
    ns.λs, ns.μs = solid.lame_parameters()
    ns.g = -solid.gravity * ns.δ[1]

    ns.d = topo.field('d', btype='std', degree=2, shape=(
        2,)) * domain.cylinder_radius  # deformation at the end of the timestep
    ns.v, ns.a = dynamic.newmark_defo_fields(ns.d)

    ns.dtest = function.replace_arguments(
        ns.d, 'd:dtest') / (solid.shear_modulus * domain.cylinder_radius**2)
    ns.traction = function.field(
        'traction', rw_sample.basis(), shape=(2,)) * precice_stress

    ns.x = ns.xref + ns.d
    ns.define_for('x', gradient='∇', jacobians=('dV', 'dS'))

    ns.F_ij = '∇ref_j(x_i)'  # deformation gradient tensor
    ns.C_ij = 'F_ki F_kj'  # right Cauchy-Green deformation tensor
    ns.E_ij = '.5 (C_ij - δ_ij)'  # Green-Lagrangian strain tensor
    ns.S_ij = 'λs E_kk δ_ij + 2 μs E_ij'  # 2nd Piola–Kirchhoff stress
    ns.P_ij = 'F_ik S_kj'  # 1st Piola–Kirchhoff stress

    res = topo.integral(
        '(∇ref_j(dtest_i) P_ij + dtest_i ρs (a_i - g_i)) dVref' @ ns, degree=4)
    res -= rw_sample.integral('dtest_i traction_i dS' @ ns)

    sqr = topo.boundary['cylinder'].integral('d_k d_k' @ ns, degree=4) / 'm2'
    cons = System(sqr, trial='d').solve_constraints(droptol=1e-10)

    system = System(res, trial='d', test='dtest')

    # initial values
    args = {a: numpy.zeros(function.arguments_for(res)[a].shape) for a in 'dt'}

    bbezier = topo.boundary.sample('bezier', 3)
    x_bbz = bbezier.bind(ns.x)

    with log.iter.plain('time step', count()) as iter_context:
        while participant.is_coupling_ongoing():

            if participant.requires_writing_checkpoint():
                checkpoint = args

            timestep = participant.get_max_time_step_size()

            # update time derivatives (implies copy)
            args = dynamic.newmark_update(timestep, **args)
            # receive boundary tractions
            args['traction'] = participant.read_data(
                rw_name, 'Stress', rw_ids, timestep)
            # solve momentum equation time step
            args = system.solve(arguments=args, constrain=cons, tol=1e-10)
            # send flap displacements
            participant.write_data(
                rw_name, 'Displacement', rw_ids, rw_sample.eval(ns.d, args) / precice_length)

            participant.advance(timestep)

            if participant.requires_reading_checkpoint():
                args = checkpoint

            if participant.is_time_window_complete():
                next(iter_context)

                t = args['t'] * precice_time

                xb = function.eval(x_bbz, args)
                with export.mplfigure('solution.jpg', dpi=150) as fig:
                    ax = fig.add_subplot(
                        111, title=f'deformation at t={t:.3s}', ylabel='[m]')
                    export.triplot(ax, xb / 'm', hull=bbezier.tri, linewidth=1)
                    ax.set_xlim(0, 4 * domain.x_center / 'm')
                    ax.set_ylim(0, 2 * domain.y_center / 'm')

                ux, uy = topo.points['A'].sample(
                    'gauss', 1).eval(ns.d, arguments=args)[0]
                log.info(f'uy: {uy:mm}')
                log.info(f'ux: {ux:mm}')
                with export.mplfigure('tip-displacement.jpg', dpi=150) as fig:
                    dynamic.add_and_plot(
                        'uy', t, 's', uy, 'mm', ax=fig.add_subplot(211))
                    dynamic.add_and_plot(
                        'ux', t, 's', ux, 'mm', ax=fig.add_subplot(212, xlabel='time [s]'))

    participant.finalize()


if __name__ == '__main__':
    cli.run(main)
