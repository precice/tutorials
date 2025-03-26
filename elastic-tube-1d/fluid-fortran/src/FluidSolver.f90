program FluidSolver
  use FluidComputeSolution, only: fluid_compute_solution
  use Utilities, only: write_vtk
  implicit none
  integer, parameter :: dp = kind(1.0d0) ! Double precision

  character(LEN=50) :: configFileName
  character(LEN=50)  :: solverName
  character(LEN=50)  :: meshName, pressureName, crossSectionLengthName
  character(LEN=50) :: outputFilePrefix
  integer             :: rank, commsize, ongoing, dimensions, bool
  integer             :: domainSize, chunkLength
  integer             :: i, j, info
  real(dp)            :: dt, t, cellwidth
  real(dp), allocatable :: pressure(:), pressure_old(:)
  real(dp), allocatable :: crossSectionLength(:), crossSectionLength_old(:)
  real(dp), allocatable :: velocity(:), velocity_old(:)
  integer, allocatable :: vertexIDs(:)
  integer             :: out_counter
  real(dp), parameter :: PI = 3.141592653589793_dp
  real(dp)            :: kappa, l
  real(dp)            :: r0, a0, u0, ampl, frequency, t_shift, p0, vel_in_0
  real(dp), allocatable :: grid(:)


  write(*, *) 'Fluid: Starting Fortran solver...'

  if (command_argument_count() /= 1) then
    write(*, *) ""
    write(*, *) "Fluid: Usage: FluidSolver <configurationFileName>"
    stop -1
  end if

  call getarg(1, configFileName)

  solverName = 'Fluid'
  outputFilePrefix = './output/out_fluid'

  ! Configure precice
  rank = 0
  commsize = 1
  call precicef_create(solverName, configFileName, rank, commsize)
  write(*, *) "preCICE configured..."

  ! Define mesh and data names
  meshName = "Fluid-Nodes-Mesh"
  pressureName = "Pressure"
  crossSectionLengthName = "CrossSectionLength"

  domainSize = 100
  chunkLength = domainSize + 1
  kappa = 100.0_dp
  l = 10.0_dp

  ! Get mesh dimensions 
  call precicef_get_mesh_dimensions(meshName, dimensions)

  ! Allocate arrays
  allocate(vertexIDs(chunkLength))
  allocate(pressure(chunkLength))
  allocate(pressure_old(chunkLength))
  allocate(crossSectionLength(chunkLength))
  allocate(crossSectionLength_old(chunkLength))
  allocate(velocity(chunkLength))
  allocate(velocity_old(chunkLength))
  allocate(grid(dimensions*chunkLength))

  ! Initialize physical parameters
  r0 = 1.0_dp / sqrt(PI) ! radius of the tube
  a0 = r0**2 * PI        ! cross-sectional area
  u0 = 10.0_dp           ! mean velocity
  ampl = 3.0_dp          ! amplitude of varying velocity
  frequency = 10.0_dp    ! frequency of variation
  t_shift = 0.0_dp       ! temporal shift of variation
  p0 = 0.0_dp            ! pressure at outlet
  vel_in_0 = u0 + ampl * sin(frequency * (t_shift) * PI)

  ! Initialize data arrays
  pressure = p0
  pressure_old = pressure
  crossSectionLength = a0
  crossSectionLength_old = crossSectionLength
  velocity = vel_in_0
  velocity_old = velocity

  ! Initialize grid coordinates
  cellwidth = l / real(domainSize, dp)
  do i = 1, chunkLength
    do j = 1, dimensions
      if (j == 1) then
        grid((i - 1)*dimensions + j) = real(i - 1, dp) * cellwidth
      else
        grid((i - 1)*dimensions + j) = 0.0_dp
      end if
    end do
  end do

  ! Initialize vertexIDs (0-based IDs)
  do i = 1, chunkLength
    vertexIDs(i) = i - 1
  end do

  call precicef_set_vertices(meshName, chunkLength, grid, vertexIDs)

  ! Check if Initial Data is Required and Write if Necessary
  call precicef_requires_initial_data(bool)
  if (bool == 1) then
    write (*, *) 'Fluid: Writing initial data'
  end if

  write (*, *) "Initialize preCICE..."
  call precicef_initialize()

  ! read initial cross-Section length
  call precicef_read_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, 0.0d0, crossSectionLength)

  ! Copy current cross-Section length to old array
  crossSectionLength_old = crossSectionLength

  ! initialize such that mass conservation is fulfilled
  do i = 1, chunkLength
    velocity_old(i) = vel_in_0*crossSectionLength_old(1)/crossSectionLength_old(i)
  end do

  t = 0.0d0 
  out_counter = 0

  ! Main coupling loop
  call precicef_is_coupling_ongoing(ongoing)
  do while (ongoing /= 0)
    ! checkpointing is required in implicit coupling
    call precicef_requires_writing_checkpoint(bool)
    if (bool .eq. 1) then
      ! nothing 
    end if

    call precicef_get_max_time_step_size(dt)

    ! solve
    call fluid_compute_solution( &
         velocity_old, pressure_old, crossSectionLength_old, &
         crossSectionLength, &
         t + dt, &          ! used for inlet velocity
         domainSize, &
         kappa, &
         dt, &              ! tau
         velocity, pressure, & ! resulting velocity pressure
         info)

    call precicef_write_data(meshName, pressureName, chunkLength, vertexIDs, pressure)
    
    call precicef_advance(dt)

    call precicef_get_max_time_step_size(dt)

    call precicef_read_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, dt, crossSectionLength)

    call precicef_requires_reading_checkpoint(bool)
    if (bool .eq. 1) then
      ! not yet converged
    else
      ! converged, advance in time 
      t = t + dt
 
      call write_vtk(t, out_counter, outputFilePrefix, chunkLength, grid, velocity, pressure, crossSectionLength)
      crossSectionLength_old = crossSectionLength
      pressure_old = pressure
      velocity_old = velocity

      out_counter = out_counter + 1
    end if

    ! Check if coupling is still ongoing
    call precicef_is_coupling_ongoing(ongoing)
  end do

  ! finalize precice and deallocate arrays 
  call precicef_finalize()
  write (*, *) 'Exiting FluidSolver'

  deallocate(pressure)
  deallocate(pressure_old)
  deallocate(crossSectionLength)
  deallocate(crossSectionLength_old)
  deallocate(velocity)
  deallocate(velocity_old)
  deallocate(grid)
  deallocate(vertexIDs)

end program FluidSolver
