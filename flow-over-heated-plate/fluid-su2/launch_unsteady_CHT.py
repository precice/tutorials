#!/usr/bin/env python


## File edited to work coupled with preCICE - Joseph Signorelli

## \file launch_unsteady_CHT_FlatPlate.py
#  \brief Python script to launch SU2_CFD with customized unsteady boundary conditions using the Python wrapper.
#  \author David Thomas
#  \version 7.4.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys
from optparse import OptionParser	# use a parser for configuration
import pysu2			            # imports the SU2 wrapped module
from math import *
import precice #import precice
import numpy
from time import sleep
# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

  # Command line options
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  # preCICE options with default settings
  parser.add_option("-p", "--precice-participant", dest="precice_name", help="Specify preCICE participant name", default="Fluid" )
  parser.add_option("-c", "--precice-config", dest="precice_config", help="Specify preCICE config file", default="../precice-config.xml")
  parser.add_option("-m", "--precice-mesh", dest="precice_mesh", help="Specify the preCICE mesh name", default="Fluid-Mesh")
  parser.add_option("-r", "--precice-read", dest="precice_read", help="Specify the preCICE read data name", default="Heat-Flux")
  parser.add_option("-w", "--precice-write", dest="precice_write", help="Specify the preCICE write data name", default="Temperature")

  (options, args) = parser.parse_args()
  options.nDim = int(2) # Specify dimension here
  options.nZone = int(1) # Specify number of zones here (1)

  # Import mpi4py for parallel run
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
  else:
    comm = 0
    rank = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
      SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return


  # Configure preCICE:
  #print("Configuring preCICE...") -- Messages appear to be hidden until sys.stdout.flush()
  size = comm.Get_size()
  try:
    interface = precice.Interface(options.precice_name, options.precice_config, rank, size)#, comm)
  except:
    print("There was an error configuring preCICE")
    return
  
  # Check preCICE + SU2 dimensions
  if options.nDim != interface.get_dimensions():
    print("SU2 and preCICE dimensions are not the same! Exiting")
    return


  CHTMarkerID = None
  CHTMarker = 'interface' # Name of CHT marker to couple

  # Get all the tags with the CHT option
  CHTMarkerList =  SU2Driver.GetAllCHTMarkersTag()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetAllBoundaryMarkers() # Returns all markers defined on this rank

  #Check if the specified marker has a CHT option and if it exists on this rank.
  if CHTMarker in CHTMarkerList and CHTMarker in allMarkerIDs.keys():
    CHTMarkerID = allMarkerIDs[CHTMarker] # So: if CHTMarkerID != None, then it exists on this rank
  
  # Number of vertices on the specified marker (per rank)
  nVertex_CHTMarker = 0         #total number of vertices (physical + halo) on this rank
  nVertex_CHTMarker_HALO = 0    #number of halo vertices
  nVertex_CHTMarker_PHYS = 0    #number of physical vertices
  iVertices_CHTMarker_PHYS = [] # indices of vertices this rank is working on
  # Datatypes must be primitive as input to SU2 wrapper code, not numpy.int8, numpy.int64, etc.. So a list is used

  # If the CHT marker is defined on this rank:
  if CHTMarkerID != None:
    nVertex_CHTMarker = SU2Driver.GetNumberVertices(CHTMarkerID) #Total number of vertices on the marker
    nVertex_CHTMarker_HALO = SU2Driver.GetNumberHaloVertices(CHTMarkerID)
    nVertex_CHTMarker_PHYS = nVertex_CHTMarker - nVertex_CHTMarker_HALO

    # Obtain indices of all vertices that are being worked on on this rank
    for iVertex in range(nVertex_CHTMarker):
      if not SU2Driver.IsAHaloNode(CHTMarkerID, iVertex):
        iVertices_CHTMarker_PHYS.append(int(iVertex))

  # Get preCICE mesh ID
  try:
    mesh_id = interface.get_mesh_id(options.precice_mesh)
  except:
    print("Invalid or no preCICE mesh name provided")
    return

  # Get coords of vertices
  coords = numpy.zeros((nVertex_CHTMarker_PHYS, options.nDim))
  for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
    coord_passive = SU2Driver.GetInitialMeshCoord(CHTMarkerID, iVertex)
    for iDim in range(options.nDim):
      coords[i, iDim] = coord_passive[iDim]

  # Set mesh vertices in preCICE:
  vertex_ids = interface.set_mesh_vertices(mesh_id, coords)

  # Get read and write data IDs
  read_data_id = interface.get_data_id(options.precice_read, mesh_id)
  write_data_id = interface.get_data_id(options.precice_write, mesh_id)

  # Instantiate arrays to hold temperature + heat flux info
  temperatures = numpy.zeros(nVertex_CHTMarker_PHYS)
  heatFluxes = numpy.zeros(nVertex_CHTMarker_PHYS)

  # Retrieve some control parameters from the driver
  deltaT = SU2Driver.GetUnsteady_TimeStep()
  TimeIter = SU2Driver.GetTime_Iter()
  nTimeIter = SU2Driver.GetnTimeIter()
  time = TimeIter*deltaT

  # Setup preCICE dt:
  precice_deltaT = interface.initialize()

  # Set up initial data for preCICE
  if (interface.is_action_required(precice.action_write_initial_data())):

    for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
      heatFluxes[i] = SU2Driver.GetVertexNormalHeatFlux(CHTMarkerID, iVertex)

    interface.write_block_scalar_data(write_data_id, vertex_ids, heatFluxes)
    interface.mark_action_fulfilled(precice.action_write_initial_data())

  interface.initialize_data()

  # Sleep briefly
  # This is critically important as I have found that initializeData is not called fast enough in CHyPS to be read here in time
  sleep(3)

  # Time loop is defined in Python so that we have access to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()


  while (interface.is_coupling_ongoing()):#TimeIter < nTimeIter):

    # Implicit coupling
    if (interface.is_action_required(precice.action_write_iteration_checkpoint())):
      # Save the state
      SU2Driver.SaveOldState()
      interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

    if (interface.is_read_data_available()):
      # Retrieve data from preCICE
      heatFluxes = interface.read_block_scalar_data(read_data_id, vertex_ids) 

      # Set the updated temperatures
      for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
          SU2Driver.SetVertexNormalHeatFlux(CHTMarkerID, iVertex, heatFluxes[i])

      # Tell the SU2 drive to update the boundary conditions
      SU2Driver.BoundaryConditionsUpdate()

    # Include barrier after preCICE stuff
    if options.with_MPI == True:
      comm.Barrier()

    # Update timestep based on preCICE
    deltaT = SU2Driver.GetUnsteady_TimeStep()
    deltaT = min(precice_deltaT, deltaT)
    SU2Driver.SetUnsteady_TimeStep(deltaT)

    # Time iteration preprocessing
    SU2Driver.Preprocess(TimeIter)

    # Run one time iteration (e.g. dual-time)
    SU2Driver.Run()

    # Postprocess the solver and exit cleanly
    SU2Driver.Postprocess()

    # Update the solver for the next time iteration
    SU2Driver.Update()
    
    # Monitor the solver and output solution to file if required
    stopCalc = SU2Driver.Monitor(TimeIter)
    
    if (interface.is_write_data_required(deltaT)):
      # Loop over the vertices
      for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
        # Get heat fluxes at each vertex
        temperatures[i] = SU2Driver.GetVertexTemperature(CHTMarkerID, iVertex)
        
      # Write data to preCICE
      interface.write_block_scalar_data(write_data_id, vertex_ids, temperatures)

    # Advance preCICE
    precice_deltaT = interface.advance(deltaT)

    # Implicit coupling:
    if (interface.is_action_required(precice.action_read_iteration_checkpoint())):
      # Reload old state
      SU2Driver.ReloadOldState()
      interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())
    else: # Output and increment as usual
      SU2Driver.Output(TimeIter)
      if (stopCalc == True):
        break
      # Update control parameters
      TimeIter += 1
      time += deltaT

    # Include barrier after preCICE stuff
    if options.with_MPI == True:
      comm.Barrier()

  interface.finalize()
  
  if SU2Driver != None:
    del SU2Driver

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
