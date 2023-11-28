! Solve material model for stresses. Written originally by Minh Hoang Nguyen
      subroutine solve_material_model(
     *   nblock, ndir, nshr, nprops,
     *   props, strains, stresses)
      ! Orientation of stresses: [11, 22, 33, 12, 23, 13]

      implicit none

      double precision :: zero
      parameter(zero=0.d0)

      integer :: i_prp_E11, i_prp_E22, i_prp_G12,
     * i_prp_nu12,i_prp_nu23

      parameter(i_prp_E11=1,
     * i_prp_E22 = 2,
     * i_prp_G12 = 3,
     * i_prp_nu12 = 4,
     * i_prp_nu23 = 5)

      integer, intent(in) :: nblock, ndir, nshr, nprops
      real*8, dimension(nprops), intent(in) :: props
      real*8, dimension(nblock, ndir + nshr), intent(in) :: strains

      real*8, dimension(nblock, ndir + nshr), intent(out) :: stresses

      real*8, dimension(ndir + nshr, ndir + nshr) :: Q
      real*8 :: E11, E22, E33, nu12, nu23, nu13, G12, G13, G23,
     * nu31, nu21, nu32
      real*8 :: temp
      integer :: k

      ! Initializations for stress calculation
      E11  = props(i_prp_E11)
      E22  = props(i_prp_E22)
      E33  = props(i_prp_E22)
      nu12 = props(i_prp_nu12)
      nu23 = props(i_prp_nu23)
      nu13 = props(i_prp_nu12)
      G12  = props(i_prp_G12)
      G13  = props(i_prp_G12)
      G23  = props(i_prp_E22)/2./(1.+props(i_prp_nu23))
      nu31 = nu13*E33/E11
      nu21 = nu12*E22/E11
      nu31 = nu13*E33/E11
      nu32 = nu23*E33/E22

      ! Calculate stresses using a material model
      temp   = nu12*nu21 + nu23*nu32 + nu13*nu31 + 2.*nu21*nu32*nu13
      Q      = zero
      Q(1,1) = (1.-nu23*nu32)*E11 / (1.-temp)
      Q(2,2) = (1.-nu13*nu31)*E22 / (1.-temp)
      Q(3,3) = (1.-nu12*nu21)*E33 / (1.-temp)
      Q(1,2) = (nu21+nu31*nu23)*E11/(1.-temp)
      Q(2,1) = Q(1,2)
      Q(1,3) = (nu31+nu21*nu32)*E11/(1.-temp)
      Q(3,1) = Q(1,3)
      Q(2,3) = (nu32+nu12*nu31)*E22/(1.-temp)
      Q(3,2) = Q(2,3)
      Q(4,4) = G12
      Q(5,5) = G23
      Q(6,6) = G13

      do k = 1, nblock
         stresses(k, :) = MATMUL(Q, strains(k, :))
      end do

      return
      end subroutine ! Subroutine solve_material_model
