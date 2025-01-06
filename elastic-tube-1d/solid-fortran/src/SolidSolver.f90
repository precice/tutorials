PROGRAM SolidSolver
  USE SolidComputeSolution_mod, ONLY: SolidComputeSolution
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(1.0D0)

  CHARACTER(LEN=512) :: configFileName
  CHARACTER(LEN=256) :: solverName
  CHARACTER(LEN=256) :: meshName, crossSectionLengthName, pressureName
  INTEGER            :: rank, commsize, ongoing, dimensions, bool
  INTEGER            :: domainSize, chunkLength
  INTEGER            :: i, j
  INTEGER, ALLOCATABLE :: vertexIDs(:)
  DOUBLE PRECISION, ALLOCATABLE :: pressure(:), crossSectionLength(:)
  DOUBLE PRECISION, ALLOCATABLE :: grid(:)
  DOUBLE PRECISION    :: dt, tubeLength, dx

  WRITE (*,*) 'Starting Solid Solver...'

  IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
    WRITE (*,*) 'Solid: Usage: SolidSolver <configurationFileName>'
    WRITE (*,*) ''
    STOP -1
  END IF

  CALL get_command_argument(1, configFileName)

  domainSize  = 4  
  chunkLength = domainSize + 1
  tubeLength  = 10.0D0

  WRITE (*,*) 'N: ', domainSize
  WRITE (*,*) 'inputs: ', COMMAND_ARGUMENT_COUNT()

  solverName              = 'Solid'
  meshName                = 'Solid-Nodes-Mesh'
  crossSectionLengthName  = 'CrossSectionLength'
  pressureName            = 'Pressure'

  rank     = 0
  commsize = 1
  CALL precicef_create(solverName, configFileName, rank, commsize)
  WRITE (*,*) 'preCICE configured...'

  CALL precicef_get_mesh_dimensions(meshName, dimensions)

  ! Allocate Arrays
  ALLOCATE(pressure(chunkLength))
  ALLOCATE(crossSectionLength(chunkLength))
  ALLOCATE(grid(dimensions * chunkLength))
  ALLOCATE(vertexIDs(chunkLength))

  pressure          = 0.0D0
  crossSectionLength = 1.0D0
  dx = tubeLength / domainSize
  DO i = 1, chunkLength
    grid((i - 1) * dimensions + 1) = dx * REAL(i - 1, DP)  ! x-coordinate
    grid((i - 1) * dimensions + 2) = 0.0D0                 ! y-coordinate
    vertexIDs(i) = i - 1                                  ! 0-based indexing here
  END DO

  ! Print the grid
  print *, "Grid values:"
  do i = 1, chunkLength
      do j = 1, dimensions
          print "(A,I4,A,F6.2)", "grid(", (i - 1) * dimensions + j, ") = ", grid((i - 1) * dimensions + j)
      end do
  end do

  CALL precicef_set_vertices(meshName, chunkLength, grid, vertexIDs)

  ! Check if Initial Data is Required and Write if Necessary
  CALL precicef_requires_initial_data(bool)
  IF (bool == 1) THEN
    CALL precicef_write_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, crossSectionLength)
  END IF

  WRITE (*,*) 'Initialize preCICE...'
  CALL precicef_initialize()

  ! Coupling Loop
  CALL precicef_is_coupling_ongoing(ongoing)
  DO WHILE (ongoing /= 0)

    CALL precicef_requires_writing_checkpoint(bool)
    IF (bool == 1) THEN
      WRITE (*,*) 'Solid: Writing iteration checkpoint (not implemented).'
    END IF

    CALL precicef_get_max_time_step_size(dt)

    CALL precicef_read_data(meshName, pressureName, chunkLength, vertexIDs, dt, pressure)

    CALL SolidComputeSolution(chunkLength, pressure, crossSectionLength)

    CALL precicef_write_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, crossSectionLength)
    
    CALL precicef_advance(dt)

    CALL precicef_requires_reading_checkpoint(bool)
    IF (bool == 1) THEN
      WRITE (*,*) 'Solid: Reading iteration checkpoint (not implemented).'
    END IF

    CALL precicef_is_coupling_ongoing(ongoing)
  END DO

  WRITE (*,*) 'Exiting SolidSolver'

  CALL precicef_finalize()

  DEALLOCATE(pressure)
  DEALLOCATE(crossSectionLength)
  DEALLOCATE(grid)
  DEALLOCATE(vertexIDs)

END PROGRAM SolidSolver
