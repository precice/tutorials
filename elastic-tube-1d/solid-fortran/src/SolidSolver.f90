program SolidSolver
  use SolidComputeSolution, only: solid_compute_solution
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  character(len=50) :: configFileName
  character(len=50) :: solverName
  character(len=50) :: meshName, crossSectionLengthName, pressureName
  integer            :: rank, commsize, ongoing, dimensions, bool
  integer            :: domainSize, chunkLength
  integer            :: i
  integer, allocatable :: vertexIDs(:)
  real(dp), allocatable :: pressure(:), crossSectionLength(:)
  real(dp), allocatable :: grid(:)
  real(dp) :: dt, tubeLength, dx

  write(*, *) 'Starting Solid Solver...'

  if (command_argument_count() /= 1) then
    write(*, *) 'Solid: Usage: SolidSolver <configurationFileName>'
    write(*, *) ''
    stop -1
  end if

  call get_command_argument(1, configFileName)

  domainSize = 100
  chunkLength = domainSize + 1
  tubeLength = 10.0_dp

  write(*, *) 'N: ', domainSize
  write(*, *) 'inputs: ', command_argument_count()

  solverName = 'Solid'
  meshName = 'Solid-Nodes-Mesh'
  crossSectionLengthName = 'CrossSectionLength'
  pressureName = 'Pressure'

  rank = 0
  commsize = 1
  call precicef_create(solverName, configFileName, rank, commsize)
  write(*, *) 'preCICE configured...'

  call precicef_get_mesh_dimensions(meshName, dimensions)

  ! Allocate arrays
  allocate(pressure(chunkLength))
  allocate(crossSectionLength(chunkLength))
  allocate(grid(dimensions*chunkLength))
  allocate(vertexIDs(chunkLength))

  pressure = 0.0_dp
  crossSectionLength = 1.0_dp
  dx = tubeLength / real(domainSize, dp)
  do i = 1, chunkLength
    grid((i - 1)*dimensions + 1) = dx * real(i - 1, dp)  ! x-coordinate
    grid((i - 1)*dimensions + 2) = 0.0_dp                 ! y-coordinate
    vertexIDs(i) = i - 1                                  ! 0-based indexing here
  end do

  call precicef_set_vertices(meshName, chunkLength, grid, vertexIDs)

  ! Check if initial data is required and write if necessary
  call precicef_requires_initial_data(bool)
  if (bool == 1) then
    call precicef_write_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, crossSectionLength)
  end if

  write (*, *) 'Initialize preCICE...'
  call precicef_initialize()

  ! Coupling loop
  call precicef_is_coupling_ongoing(ongoing)
  do while (ongoing /= 0)

    call precicef_requires_writing_checkpoint(bool)
    if (bool .eq. 1) then
      ! Do nothing here
    end if

    call precicef_get_max_time_step_size(dt)

    call precicef_read_data(meshName, pressureName, chunkLength, vertexIDs, dt, pressure)

    call solid_compute_solution(chunkLength, pressure, crossSectionLength)

    call precicef_write_data(meshName, crossSectionLengthName, chunkLength, vertexIDs, crossSectionLength)
    
    call precicef_advance(dt)

    call precicef_requires_reading_checkpoint(bool)
    if (bool .eq. 1) then
      ! nothing
    end if

    call precicef_is_coupling_ongoing(ongoing)
  end do

  write (*, *) 'Exiting SolidSolver'

  call precicef_finalize()

  deallocate(pressure)
  deallocate(crossSectionLength)
  deallocate(grid)
  deallocate(vertexIDs)

end program SolidSolver
