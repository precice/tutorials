MODULE Utilities
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
CONTAINS

  SUBROUTINE write_vtk(t, iteration, filenamePrefix, nSlices, &
                       grid, velocity, pressure, diameter)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)       :: t
    INTEGER, INTENT(IN)                :: iteration
    CHARACTER(LEN=*), INTENT(IN)       :: filenamePrefix
    INTEGER, INTENT(IN)                :: nSlices
    DOUBLE PRECISION, INTENT(IN)       :: grid(:), velocity(:), pressure(:), diameter(:)

    INTEGER :: ioUnit, i, ioStatus
    CHARACTER(LEN=256) :: filename

    WRITE(filename, '(A,"_",I0,".vtk")') TRIM(filenamePrefix), iteration
    PRINT *, 'Writing timestep at t=', t, ' to ', TRIM(filename)

    OPEN(newunit=ioUnit, file=TRIM(filename), status="replace", action="write", form="formatted", iostat=ioStatus)
    IF (ioStatus /= 0) THEN
      PRINT *, 'Error: Unable to open file ', TRIM(filename)
      RETURN
    END IF

    WRITE(ioUnit, '(A)') '# vtk DataFile Version 2.0'
    WRITE(ioUnit, '(A)') ''
    WRITE(ioUnit, '(A)') 'ASCII'
    WRITE(ioUnit, '(A)') ''
    WRITE(ioUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(ioUnit, '(A)') ''

    WRITE(ioUnit, '(A,I0,A)') 'POINTS ', nSlices, ' float'
    WRITE(ioUnit, '(A)') ''
    DO i = 1, nSlices
      WRITE(ioUnit, '(ES24.16,1X,ES24.16,1X,ES24.16)') grid(2*(i-1)+1), grid(2*(i-1)+2), 0.0D0
    END DO
    WRITE(ioUnit, '(A)') ''

    WRITE(ioUnit, '(A,I0)') 'POINT_DATA ', nSlices
    WRITE(ioUnit, '(A)') ''

    WRITE(ioUnit, '(A,A,A)') 'VECTORS ', 'velocity', ' float'
    DO i = 1, nSlices
      WRITE(ioUnit, '(ES24.16,1X,ES24.16,1X,ES24.16)') velocity(i), 0.0D0, 0.0D0
    END DO
    WRITE(ioUnit, '(A)') ''

    WRITE(ioUnit, '(A,A,A)') 'SCALARS ', 'pressure', ' float'
    WRITE(ioUnit, '(A)') 'LOOKUP_TABLE default'
    DO i = 1, nSlices
      WRITE(ioUnit, '(ES24.16)') pressure(i)
    END DO
    WRITE(ioUnit, '(A)') ''

    WRITE(ioUnit, '(A,A,A)') 'SCALARS ', 'diameter', ' float'
    WRITE(ioUnit, '(A)') 'LOOKUP_TABLE default'
    DO i = 1, nSlices
      WRITE(ioUnit, '(ES24.16)') diameter(i)
    END DO
    WRITE(ioUnit, '(A)') ''

    CLOSE(ioUnit)

  END SUBROUTINE write_vtk

END MODULE Utilities
