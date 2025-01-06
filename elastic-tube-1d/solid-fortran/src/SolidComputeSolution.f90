module SolidComputeSolution_mod
  implicit none
contains

  subroutine SolidComputeSolution(chunkLength, pressure, crossSectionLength)
    integer, intent(in) :: chunkLength
    real(8), intent(in) :: pressure(1:chunkLength)
    real(8), intent(inout) :: crossSectionLength(1:chunkLength)

    real(8) :: PI, E, r0, c_mk, c_mk2
    real(8) :: pressure0
    integer :: i

    PI = 3.14159265359d0
    E  = 10000.d0
    r0 = 1.d0 / sqrt(PI)
    c_mk = sqrt(E / (2.d0*r0))
    c_mk2 = c_mk * c_mk
    pressure0 = 0.d0

    do i=1, chunkLength
       crossSectionLength(i) = ( (pressure0 - 2.d0 * c_mk2)**2 ) / &
                               ( (pressure(i) - 2.d0 * c_mk2)**2 )
    end do

  end subroutine SolidComputeSolution

end module SolidComputeSolution_mod
