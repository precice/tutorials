program FluidSolver
  implicit none

  character(len=512) :: configFileName
  character(len=*), parameter :: solverName = 'Fluid'
  character(len=*), parameter :: meshName   = 'Fluid-Nodes-Mesh'

  integer :: rank, commsize, ongoing, bool
  double precision :: dt

  ! We “declare” the preCICE native Fortran routines as external:
  external :: precicef_create
  external :: precicef_initialize
  external :: precicef_is_coupling_ongoing
  external :: precicef_get_max_time_step_size
  external :: precicef_advance
  external :: precicef_finalize

  write(*,*) "Starting minimal Fluid Solver..."

  if (command_argument_count() /= 1) then
    write(*,*) "Usage: FluidSolver <precice-config.xml>"
    stop
  endif
  call get_command_argument(1, configFileName)

  rank = 0
  commsize = 1

  ! create participant
  call precicef_create(solverName, configFileName, rank, commsize)

  ! initialize
  call precicef_initialize()

  ! minimal coupling loop
  call precicef_is_coupling_ongoing(ongoing)
  do while(ongoing /= 0)
    call precicef_get_max_time_step_size(dt)
    call precicef_advance(dt)
    call precicef_is_coupling_ongoing(ongoing)
  end do

  ! finalize
  call precicef_finalize()

  write(*,*) "Exit FluidSolver"
end program FluidSolver
