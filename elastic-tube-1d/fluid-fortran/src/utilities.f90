module Utilities
  implicit none
  integer, parameter :: dp = kind(1.0D0)
contains

  subroutine write_vtk(t, iteration, filenamePrefix, nSlices, &
                       grid, velocity, pressure, diameter)
    implicit none

    real(dp), intent(IN)       :: t
    integer, intent(IN)        :: iteration
    character(LEN=*), intent(IN) :: filenamePrefix
    integer, intent(IN)        :: nSlices
    real(dp), intent(IN)       :: grid(:), velocity(:), pressure(:), diameter(:)

    integer :: ioUnit, i, ioStatus
    character(LEN=256) :: filename

    write (filename, '(A,"_",I0,".vtk")') trim(filenamePrefix), iteration
    print '(A, F7.6, A, A)', 'writing timestep at t=', t, ' to ', trim(filename)

    open (newunit=ioUnit, file=trim(filename), status="replace", action="write", form="formatted", iostat=ioStatus)
    if (ioStatus /= 0) then
      print *, 'Error: Unable to open file ', trim(filename)
      return
    end if

    ! Write vtk headers
    write (ioUnit, '(A)') '# vtk DataFile Version 2.0'
    write (ioUnit, '(A)') ''
    write (ioUnit, '(A)') 'ASCII'
    write (ioUnit, '(A)') ''
    write (ioUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'
    write (ioUnit, '(A)') ''

    ! Write points
    write (ioUnit, '(A,I0,A)') 'POINTS ', nSlices, ' float'
    write (ioUnit, '(A)') ''
    do i = 1, nSlices
      write (ioUnit, '(ES24.16,1X,ES24.16,1X,ES24.16)') grid(2*(i - 1) + 1), grid(2*(i - 1) + 2), 0.0D0
    end do
    write (ioUnit, '(A)') ''

    write (ioUnit, '(A,I0)') 'POINT_DATA ', nSlices
    write (ioUnit, '(A)') ''

    ! Write velocity vector field
    write (ioUnit, '(A,A,A)') 'VECTORS ', 'velocity', ' float'
    do i = 1, nSlices
      write (ioUnit, '(ES24.16,1X,ES24.16,1X,ES24.16)') velocity(i), 0.0D0, 0.0D0
    end do
    write (ioUnit, '(A)') ''

    ! Write pressure
    write (ioUnit, '(A,A,A)') 'SCALARS ', 'pressure', ' float'
    write (ioUnit, '(A)') 'LOOKUP_TABLE default'
    do i = 1, nSlices
      write (ioUnit, '(ES24.16)') pressure(i)
    end do
    write (ioUnit, '(A)') ''

    ! Write diameter
    write (ioUnit, '(A,A,A)') 'SCALARS ', 'diameter', ' float'
    write (ioUnit, '(A)') 'LOOKUP_TABLE default'
    do i = 1, nSlices
      write (ioUnit, '(ES24.16)') diameter(i)
    end do
    write (ioUnit, '(A)') ''

    close (ioUnit)

  end subroutine write_vtk

end module Utilities
