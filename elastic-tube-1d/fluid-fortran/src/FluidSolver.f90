PROGRAM FluidSolver
  USE FluidComputeSolution, ONLY: fluidComputeSolutionSerial
  USE utilities,            ONLY: write_vtk
  IMPLICIT NONE

  ! Variable Declarations
  CHARACTER(LEN=512) :: configFileName
  CHARACTER(LEN=50)  :: solverName
  CHARACTER(LEN=50)  :: meshName, pressureName, crossSectionLengthName
  CHARACTER(LEN=256) :: outputFilePrefix
  INTEGER             :: rank, commsize, ongoing, dimensions, bool
  INTEGER             :: domainSize, chunkLength
  INTEGER             :: i, j, info
  DOUBLE PRECISION    :: dt, t, cellwidth
  DOUBLE PRECISION, ALLOCATABLE :: pressure(:), pressure_old(:)
  DOUBLE PRECISION, ALLOCATABLE :: crossSectionLength(:), crossSectionLength_old(:)
  DOUBLE PRECISION, ALLOCATABLE :: velocity(:), velocity_old(:)
  INTEGER, ALLOCATABLE :: vertexIDs(:)
  INTEGER             :: out_counter
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793d0
  DOUBLE PRECISION    :: kappa, L
  DOUBLE PRECISION    :: r0, a0, u0, ampl, frequency, t_shift, p0, vel_in_0
  DOUBLE PRECISION, ALLOCATABLE :: grid(:)

  ! Start of Program
  WRITE (*,*) 'Fluid: Starting Fortran solver...'

  ! Command-Line Argument Parsing
  IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
    WRITE (*,*) ""
    WRITE (*,*) "Fluid: Usage: FluidSolver <configurationFileName>"
    STOP -1
  END IF

  CALL getarg(1, configFileName)

  solverName = 'Fluid'
  outputFilePrefix = './output/out_fluid'

  ! Initialize preCICE Interface
  rank = 0
  commsize = 1
  CALL precicef_create(solverName, configFileName, rank, commsize)
  WRITE (*,*) "preCICE configured..."

  ! Define Mesh and Data Names
  meshName = "Fluid-Nodes-Mesh"
  pressureName = "Pressure"
  crossSectionLengthName = "CrossSectionLength"

  domainSize = 100
  chunkLength = domainSize + 1
  kappa = 100.0d0
  L = 10.0d0

  ! Get Mesh Dimensions from preCICE
  CALL precicef_get_mesh_dimensions(meshName, dimensions)

  ! Allocate Arrays
  ALLOCATE(vertexIDs(chunkLength))
  ALLOCATE(pressure(chunkLength))
  ALLOCATE(pressure_old(chunkLength))
  ALLOCATE(crossSectionLength(chunkLength))
  ALLOCATE(crossSectionLength_old(chunkLength))
  ALLOCATE(velocity(chunkLength))
  ALLOCATE(velocity_old(chunkLength))
  ALLOCATE(grid(dimensions * chunkLength))

  ! Initialize vertexIDs (0-based IDs)
  DO i = 1, chunkLength
    vertexIDs(i) = i - 1
  END DO

  ! Initialize Physical Parameters
  r0 = 1.0d0 / SQRT(PI)
  a0 = r0**2 * PI
  u0 = 10.0d0
  ampl = 3.0d0
  frequency = 10.0d0
  t_shift = 0.0d0
  p0 = 0.0d0
  vel_in_0 = u0 + ampl * SIN(frequency * (t_shift) * PI)

  ! Initialize Data Arrays
  pressure = p0
  pressure_old = pressure
  crossSectionLength = a0
  crossSectionLength_old = crossSectionLength
  velocity = vel_in_0
  velocity_old = velocity

  ! Initialize Grid Coordinates
  cellwidth = L / REAL(domainSize, KIND=8)
  DO i = 1, chunkLength
    DO j = 1, dimensions
      IF (j == 1) THEN
        grid((i - 1) * dimensions + j) = REAL(i - 1, KIND=8) * cellwidth
      ELSE
        grid((i - 1) * dimensions + j) = 0.0d0
      END IF
    END DO
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
    WRITE (*,*) 'Fluid: Writing initial data'
  END IF

  ! Initialize Simulation Time
  t = 0.0d0
  WRITE (*,*) "Initialize preCICE..."
  CALL precicef_initialize()

  ! Read Initial Cross-Section Length
  CALL precicef_read_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, 0.0d0, crossSectionLength)

  ! Copy Current Cross-Section Length to Old Array
  crossSectionLength_old = crossSectionLength

  ! initialize such that mass conservation is fulfilled
  DO i = 1, chunkLength
    velocity_old(i) = vel_in_0 * crossSectionLength_old(1) / crossSectionLength_old(i)
  END DO

  ! Initialize Output Counter
  out_counter = 0

  ! Print all arrays with 2 decimal places
  print *, "Index | Pressure | Pressure_Old | CrossSection | CrossSection_Old | Velocity | Velocity_Old"
  do i = 1, chunkLength
      print "(I5, 3X, F8.2, 3X, F13.2, 3X, F13.2, 3X, F16.2, 3X, F8.2, 3X, F13.2)", &
          i - 1, pressure(i), pressure_old(i), crossSectionLength(i), &
          crossSectionLength_old(i), velocity(i), velocity_old(i)
  end do 

  ! Main Coupling Loop
  CALL precicef_is_coupling_ongoing(ongoing)
  DO WHILE (ongoing /= 0)
    ! Check if Writing a Checkpoint is Required
    CALL precicef_requires_writing_checkpoint(bool)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'Fluid: Writing iteration checkpoint'
    END IF

    ! Get Maximum Time Step Size from preCICE
    CALL precicef_get_max_time_step_size(dt)

    ! Compute Fluid Solution
    CALL fluidComputeSolutionSerial( &
         velocity_old, pressure_old, crossSectionLength_old, &
         crossSectionLength, &
         t + dt, &          ! used for inlet velocity
         domainSize, &    
         kappa, &        
         dt, &              ! tau
         velocity, pressure, & ! resulting velocity pressure
         info) 

    ! Print all arrays with 2 decimal places
    print *, "Index | Pressure | Pressure_Old | CrossSection | CrossSection_Old | Velocity | Velocity_Old"
    do i = 1, chunkLength
        print "(I5, 3X, F8.2, 3X, F13.2, 3X, F13.2, 3X, F16.2, 3X, F8.2, 3X, F13.2)", &
            i - 1, pressure(i), pressure_old(i), crossSectionLength(i), &
            crossSectionLength_old(i), velocity(i), velocity_old(i)
    end do 

    CALL precicef_write_data(meshName, pressureName, chunkLength, vertexIDs, pressure)
    
    CALL precicef_advance(dt)

    CALL precicef_get_max_time_step_size(dt)

    CALL precicef_read_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, dt, crossSectionLength)

    CALL precicef_requires_reading_checkpoint(bool)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'Fluid: Reading iteration checkpoint'
    ELSE
      t = t + dt
 
      CALL write_vtk(t, out_counter, outputFilePrefix, chunkLength, grid, velocity, pressure, crossSectionLength)
      crossSectionLength_old = crossSectionLength
      pressure_old = pressure
      velocity_old = velocity

      out_counter = out_counter + 1
    END IF

    ! Check if Coupling is Still Ongoing
    CALL precicef_is_coupling_ongoing(ongoing)
  END DO

  ! Finalize preCICE Interface
  CALL precicef_finalize()
  WRITE (*,*) 'Exiting FluidSolver'

  ! Deallocate Dynamically Allocated Arrays
  DEALLOCATE(pressure)
  DEALLOCATE(pressure_old)
  DEALLOCATE(crossSectionLength)
  DEALLOCATE(crossSectionLength_old)
  DEALLOCATE(velocity)
  DEALLOCATE(velocity_old)
  DEALLOCATE(grid)
  DEALLOCATE(vertexIDs)

END PROGRAM FluidSolver
