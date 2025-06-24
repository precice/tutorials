#!/usr/bin/env python3

import numpy as np
import precice
import os
from scipy.interpolate import interp1d

import ufl
from dune.fem.space import lagrange as solutionSpace
from dune.ufl import DirichletBC, Constant
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.operator import galerkin
from dune.fem.utility import Sampler
from dune.grid import cartesianDomain
from dune.alugrid import aluGrid
from dune.fem.function import uflFunction, gridFunction
from dune.ufl import expression2GF

if __name__ == '__main__':

    # standard heat equation: u_t + k \Delta u = 0
    k = 100  # kg * m / s^3 / K, https://en.wikipedia.org/wiki/Thermal_conductivity

    # Create mesh and define function space
    nx = 100
    ny = 25

    dt_out = 0.2  # interval for writing VTK files
    y_top = 0
    y_bottom = y_top - .25
    x_left = 0
    x_right = x_left + 1

    # preCICE setup
    participant = precice.Participant("Solid", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Solid-Mesh"
    interface_x_coordinates = np.linspace(
        0, 1, nx + 1)  # meshpoints for interface values
    vertices = [[x0, 0] for x0 in interface_x_coordinates]

    vertex_ids = participant.set_mesh_vertices(mesh_name, vertices)

    domain = cartesianDomain([x_left, y_bottom], [x_right, y_top], [nx, ny])
    mesh = aluGrid(domain, 2, 2, elementType="simplex")  # create a simple mesh with dimGrid=2 and dimWorld=2
    space = solutionSpace(mesh, order=1)
    u = ufl.TrialFunction(space)
    v = ufl.TestFunction(space)
    x = ufl.SpatialCoordinate(ufl.triangle)
    n = ufl.FacetNormal(space)

    u0 = uflFunction(mesh, name="u0", order=space.order,
                     ufl=Constant(310))
    uold = space.interpolate(u0, name='uold')
    unew = space.interpolate(u0, name='unew')
    ucheckpoint = uold.copy(name='u_checkpoint')

    # creating boundary condition with data from the fluid

    ug_interp = interp1d(interface_x_coordinates, np.zeros(nx + 1))
    def ug_f(x): return ug_interp(x[0])

    @gridFunction(mesh, name="u_gamma", order=1)
    def u_gamma(xg):
        return ug_f(xg)

    eps = 1e-8
    bc_top = DirichletBC(space, u_gamma, x[1] > - eps)
    bc_bottom = DirichletBC(space, Constant(310.), x[1] < -0.25 + eps)
    bcs = [bc_bottom, bc_top]

    k = Constant(k)
    dt = Constant(0.01, name='dt')

    A = u * v * ufl.dx + dt * k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    b = uold * v * ufl.dx

    scheme = solutionScheme([A == b, *bcs], solver='cg')

    # Weak form of the flux
    flux_expr = -(u - uold) * v / dt * ufl.dx - k * \
        ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    flux_expr_operator = galerkin(flux_expr)
    flux_sol = space.interpolate(u0, name='flux_sol')

    # Sample the flux on the top edge
    flux_sol_expr = expression2GF(
        flux_sol.space.grid, flux_sol, flux_sol.space.order)
    sampler_weak_flux = Sampler(flux_sol_expr)

    def flux_f_weak(): return sampler_weak_flux.lineSample(
        [0., 0.], [1., 0.], nx + 1)[1] * nx

    if not os.path.exists("output"):
        os.makedirs("output")
    vtk = mesh.sequencedVTK("output/heat", pointdata=[unew])

    participant.initialize()
    precice_dt = participant.get_max_time_step_size()
    dt.assign(min(float(dt), precice_dt))

    t = float(dt)

    while participant.is_coupling_ongoing():

        # write checkpoint
        if participant.requires_writing_checkpoint():
            t_check = t
            ucheckpoint.assign(uold)
        ug = participant.read_data(mesh_name, "Temperature", vertex_ids, dt)
        ug_interp.y = ug

        scheme.model.dt = dt
        scheme.solve(target=unew)

        # compute flux, solution goes into flux_sol
        flux_expr_operator(unew, flux_sol)
        flux_values = flux_f_weak()  # sample the flux function

        participant.write_data(mesh_name, "Heat-Flux", vertex_ids, flux_values)

        participant.advance(dt)
        precice_dt = participant.get_max_time_step_size()
        dt.assign(min(float(dt), precice_dt))

        # roll back to checkpoint
        if participant.requires_reading_checkpoint():
            t = t_check
            uold.assign(ucheckpoint)

        else:  # update solution

            uold.assign(unew)
            t += float(dt)

        if participant.is_time_window_complete():
            # we need some tolerance, since otherwise output might be skipped.
            tol = 10e-5
            if abs((t + tol) % dt_out) < 2 * tol:  # output if t is a multiple of dt_out
                print("output vtk for time = {}".format(float(t)))
                vtk()

    participant.finalize()
