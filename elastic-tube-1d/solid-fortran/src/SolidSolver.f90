program SolidSolver
  implicit none

  character(len=512) :: configFileName
  character(len=*), parameter :: solverName = 'Solid'

  integer :: rank, commsize, ongoing
  double precision :: dt

  external :: precicef_create
  external :: precicef_initialize
  external :: precicef_is_coupling_ongoing
  external :: precicef_get_max_time_step_size
  external :: precicef_advance
  external :: precicef_finalize

  write(*,*) "Starting minimal Solid Solver..."

  if(command_argument_count() /= 1) then
    write(*,*) "Usage: SolidSolver <precice-config.xml>"
    stop
  end if
  call get_command_argument(1, configFileName)

  rank = 0
  commsize = 1

  call precicef_create(solverName, configFileName, rank, commsize)
  
  call precicef_initialize()

  call precicef_is_coupling_ongoing(ongoing)
  do while(ongoing /= 0)
    call precicef_get_max_time_step_size(dt)
    call precicef_advance(dt)
    call precicef_is_coupling_ongoing(ongoing)
  end do

  call precicef_finalize()

  write(*,*) "Exit SolidSolver"
end program SolidSolver
