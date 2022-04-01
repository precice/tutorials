#!/usr/bin/env python3

import numpy as np
import precice

from scipy.interpolate import interp1d

import ufl
from dune.fem.space import lagrange as solutionSpace
from dune.ufl import DirichletBC, Constant
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.operator import galerkin
from dune.fem.utility import Sampler, lineSample 
from dune.grid import cartesianDomain
from dune.grid import structuredGrid as leafGridView
from dune.alugrid import aluSimplexGrid
from dune.fem.function import uflFunction, gridFunction
from dune.ufl import expression2GF

## standard heat equation: u_t + k \Delta u = 0
class Problem_heat:
    eps = 1e-8
    ## gridsize = number of internal variables
    def __init__(self, gridsize = 20, k = 100, xa = 0., xb = 1., ya = -0.25, yb = 0.0):
        xa, xb, ya, yb = float(xa), float(xb), float(ya), float(yb)
        ## n = gridsize = number of internal unknowns
        ## => n + 2 nodes per unit length
        ## => n + 1 cells per unit length
        self.NN = gridsize + 2
        self.xx = np.linspace(0, 1, self.NN) ## meshpoints for interface values
        self.dx = 1./(gridsize + 1)

        self.k = Constant(k)
        self.dt = Constant(0., name = 'dt')

        self.domain = cartesianDomain([xa, ya], [xb, yb], [gridsize + 1, gridsize + 1])

        self.mesh = aluSimplexGrid(self.domain, serial = True)
        self.space = solutionSpace(self.mesh, order = 1, storage = 'petsc')

        self.x = ufl.SpatialCoordinate(ufl.triangle)
        
        self.u0 = uflFunction(self.mesh, name = "u0", order = self.space.order,
                              ufl = Constant(310))
        
        self.u = ufl.TrialFunction(self.space)
        self.v = ufl.TestFunction(self.space)
        self.uold = self.space.interpolate(self.u0, name = 'uold')
        self.unew = self.space.interpolate(self.u0, name = 'unew')
        self.u_checkpoint = self.uold.copy(name = 'u_checkpoint')
        
        self.A = self.u*self.v*ufl.dx + self.dt*self.k*ufl.dot(ufl.grad(self.u), ufl.grad(self.v))*ufl.dx
        self.b = self.uold*self.v*ufl.dx
        
        
        
        self.ug_interp = interp1d(self.xx, np.zeros(self.NN))

        self.ug_f = lambda x : self.ug_interp(x[0])
        
        @gridFunction(self.mesh, name="u_gamma", order=1)
        def u_gamma(xg):
            return self.ug_f(xg)
        
        self.bc_bottom = DirichletBC(self.space, Constant(310.), self.x[1] < -0.25 + self.eps)
        self.bc_top = DirichletBC(self.space, u_gamma, self.x[1] > - self.eps)
        self.bcs = [self.bc_bottom, self.bc_top]

        self.scheme = solutionScheme([self.A == self.b, *self.bcs], solver = 'cg')
        
        self.flux_expr = ((self.u - self.uold)*self.v/self.dt*ufl.dx + 
                          self.k*ufl.dot(ufl.grad(self.u), ufl.grad(self.v))*ufl.dx)
        self.flux_expr_operator = galerkin(self.flux_expr)
        self.flux_sol = self.space.interpolate(self.u0, name = 'flux_sol')
        
        self.flux_sol_expr = expression2GF(self.flux_sol.space.grid, self.flux_sol, self.flux_sol.space.order)
        self.sampler_weak_flux = Sampler(self.flux_sol_expr)
        #self.flux_f_weak = lambda : self.sampler_weak_flux.lineSample([0. + self.eps, 0. - self.eps ], [1. - self.eps, 0. - self.eps], self.NN)[1]
        self.flux_f_weak = lambda : self.sampler_weak_flux.lineSample([0., 0. - 1e-6], [1., 0. - 1e-6], self.NN)[1]
               
            
        #self.create_checkpoint()

       
    def get_flux(self):
        ## compute flux, solution goes into self.flux_sol
        self.flux_expr_operator(self.unew, self.flux_sol)
        flux = self.flux_f_weak() / self.dx 
        #flux = lineSample(self.flux_sol_expr, [0., 0.], [1., 0.], gridsize+1)[1] / self.dx
        flux[0], flux[-1] = 0, 0
        return flux        

    def load_checkpoint(self):
        self.uold.interpolate(self.u_checkpoint)
        self.unew.interpolate(self.u_checkpoint)

    def create_checkpoint(self):
        self.u_checkpoint.interpolate(self.uold)

    def do_step(self, dt, ug):
        self.scheme.model.dt = dt

        self.ug_interp.y = ug

        self.scheme.solve(target = self.unew)
        flux = self.get_flux()
        self.uold.assign(self.unew)
        return flux


if __name__ == '__main__':


    solver = Problem_heat()
    vtk = solver.mesh.sequencedVTK("heat", pointdata=[solver.unew])
    
    # preCICE setup
    interface = precice.Interface("Solid", "../precice-config.xml", 0, 1)

    # define coupling mesh
    mesh_name = "Solid-Mesh"
    mesh_id = interface.get_mesh_id(mesh_name)
    vertices = [[x0, 0] for x0 in solver.xx]
    print(len(vertices))
    vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)
    temperature_id = interface.get_data_id("Temperature", mesh_id)
    flux_id = interface.get_data_id("Heat-Flux", mesh_id)
    
    precice_dt = interface.initialize()
    
    dt = 0.01
    
    for i in range(10):
        print("Step " + str(i))
        vtk()
        
        temperature_values = interface.read_block_scalar_data(temperature_id, vertex_ids)
        
        flux_values = solver.do_step(dt, temperature_values)
        
        print(flux_values.shape)
        
        print(flux_values)
        
        interface.write_block_scalar_data(flux_id, vertex_ids, flux_values)
        
        dt = min(dt, precice_dt)
        
        precice_dt = interface.advance(dt)
        
    interface.finalize()
    
    print("Done")

    
