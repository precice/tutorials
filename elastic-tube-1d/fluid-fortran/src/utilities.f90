module utilities
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains

  subroutine write_vtk(t, iteration, filenamePrefix, nSlices, &
                       grid, velocity, pressure, diameter)
    implicit none

    ! Input Parameters
    real(dp), intent(in) :: t
    integer, intent(in) :: iteration
    character(len=*), intent(in) :: filenamePrefix
    integer, intent(in) :: nSlices
    real(dp), intent(in) :: grid(:)
    real(dp), intent(in) :: velocity(:)
    real(dp), intent(in) :: pressure(:)
    real(dp), intent(in) :: diameter(:)

    ! Local Variables
    integer :: ioUnit, i, ioStatus
    character(len=256) :: filename

    ! Generate Filename
    write(filename, '(A,"_",I0,".vtk")') trim(filenamePrefix), iteration

    ! Print Status Message
    print *, 'Writing timestep at t=', t, ' to ', trim(filename)

    ! Open File
    open(newunit=ioUnit, file=trim(filename), status="replace", &
         action="write", form="formatted", iostat=ioStatus)
    if (ioStatus /= 0) then
      print *, 'Error: Unable to open file ', trim(filename)
      return
    end if

    ! Write VTK Header
    write(ioUnit, '(A)') '# vtk DataFile Version 2.0'
    write(ioUnit, '(A)') ''
    write(ioUnit, '(A)') 'ASCII'
    write(ioUnit, '(A)') ''
    write(ioUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'
    write(ioUnit, '(A)') ''

    ! Write POINTS Section
    write(ioUnit, '(A,I0,A)') 'POINTS ', nSlices, ' float'
    write(ioUnit, '(A)') ''

    ! Export Mesh Points
    do i = 1, nSlices
      write(ioUnit, '(F24.16,1X,F24.16,1X,F24.16)') grid(2*(i-1)+1), grid(2*(i-1)+2), 0.0_dp
    end do
    write(ioUnit, '(A)') ''

    ! Write POINT_DATA Section
    write(ioUnit, '(A,I0)') 'POINT_DATA ', nSlices
    write(ioUnit, '(A)') ''

    ! Export Velocity Vectors
    write(ioUnit, '(A,A,A)') 'VECTORS ', 'velocity', ' float'
    do i = 1, nSlices
      write(ioUnit, '(F24.16,1X,F24.16,1X,F24.16)') velocity(i), 0.0_dp, 0.0_dp
    end do
    write(ioUnit, '(A)') ''

    ! Export Pressure Scalars
    write(ioUnit, '(A,A,A)') 'SCALARS ', 'pressure', ' float'
    write(ioUnit, '(A)') 'LOOKUP_TABLE default'
    do i = 1, nSlices
      write(ioUnit, '(F24.16)') pressure(i)
    end do
    write(ioUnit, '(A)') ''

    ! Export Diameter Scalars
    write(ioUnit, '(A,A,A)') 'SCALARS ', 'diameter', ' float'
    write(ioUnit, '(A)') 'LOOKUP_TABLE default'
    do i = 1, nSlices
      write(ioUnit, '(F24.16)') diameter(i)
    end do
    write(ioUnit, '(A)') ''

    ! Close File
    close(ioUnit)

  end subroutine write_vtk

end module utilities
