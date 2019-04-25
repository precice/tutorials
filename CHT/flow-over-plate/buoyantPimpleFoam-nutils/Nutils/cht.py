#! /usr/bin/env python3

import nutils, numpy, treelog
import precice
from mpi4py import MPI


def main(elemsize: 'mesh width in x and y direction' = 0.05,
            btype: 'type of basis function (std/spline)' = 'std',
            degree: 'polynomial degree' = 1,
            dt = .01):

  print("Running utils")
  
  # the mesh
  grid = [numpy.linspace(a, b, round((b-a)/size)+1) for (a,b,size) in [(0,1,elemsize), (-.25,0,elemsize), (0,.05,.05)]]  
  domain, geom = nutils.mesh.rectilinear(grid, periodic=[2])

  # nutils namespace
  ns = nutils.function.Namespace()
  ns.x = geom
  ns.basis = domain.basis(btype, degree=degree)
  ns.u = 'basis_n ?lhs_n' # solution 
  ns.dudt = 'basis_n (?lhs_n - ?lhs0_n) / ?dt'
  ns.flux = 'basis_n ?fluxdofs_n'
  ns.k = 100 # thermal diffusivity
  ns.uwall = 310 # wall temperature

  # the weak form
  res = domain.integral('(basis_n dudt + k basis_n,i u_,i) d:x' @ ns, degree=degree*2)

  # Dirichlet boundary condition
  sqr = domain.boundary['bottom'].integral('(u - uwall)^2 d:x' @ ns, degree=degree*2)
  cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15)
  
  # preCICE setup
  configFileName = "../precice-config.xml"
  participantName = "Nutils"
  solverProcessIndex = 0
  solverProcessSize = 1
  interface = precice.Interface(participantName, solverProcessIndex, solverProcessSize)
  interface.configure(configFileName)

  # define coupling meshes
  meshNameGP = "Nutils-Mesh-GP" # Gauss points
  meshNameCC = "Nutils-Mesh-CC" # cell centers (potentially sub-sampled)
  meshIDGP = interface.get_mesh_id(meshNameGP)
  meshIDCC = interface.get_mesh_id(meshNameCC)

  couplinginterface = domain.boundary['top']
  couplingsampleGP = couplinginterface.sample('gauss', degree=degree*2)
  couplingsampleCC = couplinginterface.sample('uniform', 4) # number of sub-samples for better mapping

  verticesGP = couplingsampleGP.eval(ns.x).ravel()
  verticesCC = couplingsampleCC.eval(ns.x).ravel()
  dataIndicesGP = numpy.zeros(couplingsampleGP.npoints)
  dataIndicesCC = numpy.zeros(couplingsampleCC.npoints)
  interface.set_mesh_vertices(meshIDGP, couplingsampleGP.npoints, verticesGP, dataIndicesGP)
  interface.set_mesh_vertices(meshIDCC, couplingsampleCC.npoints, verticesCC, dataIndicesCC)

  # coupling data
  writeData = "Heat-Flux"
  readData =  "Temperature"
  writedataID = interface.get_data_id(writeData, meshIDCC)
  readdataID = interface.get_data_id(readData, meshIDGP)

  # heat flux computation
  projectionmatrix = couplinginterface.integrate(ns.eval_nm('basis_n basis_m d:x'), degree=degree*2)
  projectioncons = numpy.zeros(res.shape)
  projectioncons[projectionmatrix.rowsupp(1e-15)] = numpy.nan
  fluxdofs = lambda v: projectionmatrix.solve(v, constrain=projectioncons)

  precice_dt = interface.initialize()

  cons0 = cons # to not lose the Dirichlet BC at the bottom 
  lhs0 = numpy.zeros(res.shape)
  timestep = 1

  # project initial condition and visualize
  sqr = domain.integral('(u - uwall)^2' @ ns, degree=degree*2)
  lhs0 = nutils.solver.optimize('lhs', sqr)
  bezier = domain.sample('bezier', 2)
  x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs0)
  with treelog.add(treelog.DataLog()):
    nutils.export.vtk('Solid_0', bezier.tri, x, T=u)
  

  while interface.is_coupling_ongoing():
  
    # read temperature from interface
    if interface.is_read_data_available():  
      readdata = numpy.zeros(couplingsampleGP.npoints)
      interface.read_block_scalar_data(readdataID, couplingsampleGP.npoints, dataIndicesGP, readdata)
      coupledata = couplingsampleGP.asfunction(readdata)

      sqr = couplingsampleGP.integral((ns.u - coupledata)**2)
      cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons0)
    
    # save checkpoint
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
      lhscheckpoint = lhs0
      interface.fulfilled_action(precice.action_write_iteration_checkpoint())
      
    # potentially adjust non-matching timestep sizes  
    dt = min(dt, precice_dt)  
    
    # solve nutils timestep
    lhs = nutils.solver.solve_linear('lhs', res, constrain=cons, arguments=dict(lhs0=lhs0, dt=dt))

    # write heat fluxes to interface
    if interface.is_write_data_required(dt):
      fluxvalues = res.eval(lhs0=lhs0, lhs=lhs, dt=dt)
      writedata = couplingsampleCC.eval('-flux' @ ns, fluxdofs=fluxdofs(fluxvalues))
      interface.write_block_scalar_data(writedataID, couplingsampleCC.npoints, dataIndicesCC, writedata)

    # do the coupling
    precice_dt = interface.advance(dt)

    # read checkpoint if required
    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
      interface.fulfilled_action(precice.action_read_iteration_checkpoint())
      lhs0 = lhscheckpoint
    else: # go to next timestep and visualize
      bezier = domain.sample('bezier', 2)
      x, u = bezier.eval(['x_i', 'u'] @ ns, lhs=lhs0)
      with treelog.add(treelog.DataLog()):
        if timestep % 20 == 0:
          nutils.export.vtk('Solid_' + str(timestep), bezier.tri, x, T=u)
      timestep += 1
      lhs0 = lhs


  interface.finalize()
  
if __name__ == '__main__':
  nutils.cli.run(main)
  
  
  
  
