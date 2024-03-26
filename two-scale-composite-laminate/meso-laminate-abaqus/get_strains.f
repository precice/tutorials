! Get strains from Abaqus
      subroutine get_strains(
     *   nblock, ndir, nshr, nstatev,
     *   strainInc, stateOld, state, strains)
         ! Orientation of strains: [11, 22, 33, 12, 23, 13]
 
      implicit none

      integer, intent(in) :: nblock, ndir, nshr, nstatev
      real*8, dimension(nblock, nstatev), intent(in) :: stateOld
      real*8, dimension(nblock, ndir + nshr), intent(in) :: strainInc

      real*8, dimension(nblock, nstatev), intent(out) :: state
      real*8, dimension(nblock, ndir + nshr), intent(out) :: strains

      real*8, dimension(ndir + nshr) :: strain
      integer :: k, d
      integer :: i_sdv_eps11, i_sdv_eps22, i_sdv_eps33,
     *   i_sdv_gamma12, i_sdv_gamma23, i_sdv_gamma13,
     *   i_sdv_active, i_sdv_t_f

      parameter(i_sdv_eps11=1,
     *    i_sdv_eps22 = 2,
     *    i_sdv_eps33 = 3,
     *    i_sdv_gamma12 = 4,
     *    i_sdv_gamma23 = 5,
     *    i_sdv_gamma13 = 6,
     *    i_sdv_active = 7,
     *    i_sdv_t_f = 8)

      double precision :: zero, half, one, two, three

      parameter(zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, three=3.d0)

      state(:, i_sdv_eps11) = zero
      state(:, i_sdv_eps22) = zero
      state(:, i_sdv_eps33) = zero
      state(:, i_sdv_gamma12) = zero
      state(:, i_sdv_gamma23) = zero
      state(:, i_sdv_gamma13) = zero
      state(:, i_sdv_active) = one
      state(:, i_sdv_t_f) = zero
 
      ! Loop through material points to collect strains
      do k = 1, nblock
         state(k, :) = stateOld(k, :)
 
         strain(1) = state(k, i_sdv_eps11) + strainInc(k, 1)
         strain(2) = state(k, i_sdv_eps22) + strainInc(k, 2)
         strain(3) = state(k, i_sdv_eps33) + strainInc(k, 3)
         strain(4) = state(k, i_sdv_gamma12) + two*strainInc(k, 4)
         strain(5) = state(k, i_sdv_gamma23) + two*strainInc(k, 5)
         strain(6) = state(k, i_sdv_gamma13) + two*strainInc(k, 6)
         state(k, i_sdv_eps11:i_sdv_gamma13) = strain
 
         do d = 1, ndir + nshr

            strains(k, d) = strain(d)
 
         end do ! ndir + nshr
 
      end do ! nblock
 
      return
      end subroutine ! subroutine get_strains
