module FluidComputeSolution
    implicit none
    integer, parameter :: dp = kind(1.0d0)

contains

    subroutine fluid_compute_solution(velocity_old, pressure_old, &
        crossSectionLength_old, crossSectionLength, t, N, kappa, tau, &
        velocity, pressure, info)

        real(dp), intent(in) :: velocity_old(:), pressure_old(:)
        real(dp), intent(in) :: crossSectionLength_old(:), crossSectionLength(:)
        real(dp), intent(in) :: t
        integer, intent(in) :: N
        real(dp), intent(in) :: kappa, tau
        real(dp), intent(inout) :: velocity(:), pressure(:)
        integer, intent(out) :: info

        ! Local variables
        integer :: i, k
        real(dp), parameter :: PI = 3.141592653589793_dp
        real(dp), parameter :: e = 10000.0_dp
        real(dp), parameter :: c_mk2 = e / 2.0_dp * sqrt(PI)
        real(dp), parameter :: u0 = 10.0_dp, ampl = 3.0_dp, frequency = 10.0_dp, &
                                t_shift = 0.0_dp
        real(dp), parameter :: tolerance = 1.0e-15_dp
        integer, parameter :: max_iterations = 50

        real(dp) :: alpha, L, dx, velocity_in, tmp2, norm_1, norm_2, norm

        ! LAPACK Variables
        integer :: nlhs, nrhs
        real(dp), allocatable :: Res(:)
        real(dp), allocatable :: LHS(:, :)
        integer, allocatable :: ipiv(:)

        nlhs = 2*N + 2
        nrhs = 1

        ! Allocate arrays
        allocate (Res(2*N + 2))
        allocate (LHS(2*N + 2, 2*N + 2))
        allocate (ipiv(nlhs))

        velocity = velocity_old
        pressure = pressure_old

        ! Stabilization intensity
        alpha = 0.0  !(N * kappa * tau) / (N * tau + 1);
        L = 10.0
        dx = L/kappa  !1.0 / (N * kappa);

        ! Output status from dgesv (0 = success, < 0 = invalid argument, > 0 = singular matrix)
        info = 0

        ! Nonlinear solver loop
        do k = 1, max_iterations
            ! Initialize residual vector
            Res = 0.0

            ! Compute residuals
            do i = 2, N  ! Adjusted for 1-based indexing
                ! Momentum
                Res(i) = (velocity_old(i)*crossSectionLength_old(i) - velocity(i)*crossSectionLength(i))*dx/tau
                Res(i) = Res(i) + 0.25*(-crossSectionLength(i + 1)*velocity(i)*velocity(i + 1) - &
                                        crossSectionLength(i)*velocity(i)*velocity(i + 1))
                Res(i) = Res(i) + 0.25*(-crossSectionLength(i + 1)*velocity(i)**2 - &
                                        crossSectionLength(i)*velocity(i)**2 + &
                                        crossSectionLength(i)*velocity(i - 1)*velocity(i) + &
                                        crossSectionLength(i - 1)*velocity(i - 1)*velocity(i))
                Res(i) = Res(i) + 0.25*(crossSectionLength(i - 1)*velocity(i - 1)**2 + &
                                        crossSectionLength(i)*velocity(i - 1)**2)
                Res(i) = Res(i) + 0.25*(crossSectionLength(i - 1)*pressure(i - 1) + &
                                        crossSectionLength(i)*pressure(i - 1) - &
                                        crossSectionLength(i - 1)*pressure(i) + &
                                        crossSectionLength(i + 1)*pressure(i) - &
                                        crossSectionLength(i)*pressure(i + 1) - &
                                        crossSectionLength(i + 1)*pressure(i + 1))

                ! Continuity
                Res(i + N + 1) = (crossSectionLength_old(i) - crossSectionLength(i))*dx/tau
                Res(i + N + 1) = Res(i + N + 1) + 0.25*(crossSectionLength(i - 1)*velocity(i - 1) + &
                                                crossSectionLength(i)*velocity(i - 1) + &
                                                crossSectionLength(i - 1)*velocity(i) - &
                                                crossSectionLength(i + 1)*velocity(i) - &
                                                crossSectionLength(i)*velocity(i + 1) - &
                                                crossSectionLength(i + 1)*velocity(i + 1))
                Res(i + N + 1) = Res(i + N + 1) + alpha*(pressure(i - 1) - 2.0*pressure(i) + pressure(i + 1))
            end do

            ! Boundary conditions
            velocity_in = u0 + ampl*sin(frequency*(t + t_shift)*PI)
            Res(1) = velocity_in - velocity(1)
            ! Pressure Inlet is linearly interpolated
            Res(N + 2) = -pressure(1) + 2.0*pressure(2) - pressure(3)
            ! Velocity Outlet is linearly interpolated
            Res(N + 1) = -velocity(N + 1) + 2.0*velocity(N) - velocity(N - 1)
            ! Pressure Outlet is "non-reflecting"
            tmp2 = sqrt(c_mk2 - pressure_old(N + 1)/2.0) - &
                (velocity(N + 1) - velocity_old(N + 1))/4.0
            Res(2*N + 2) = -pressure(N + 1) + 2.0*(c_mk2 - tmp2**2)

            ! Compute residual norm
            norm_1 = sqrt(sum(Res**2))
            norm_2 = sqrt(sum(pressure**2) + sum(velocity**2))
            norm = norm_1/norm_2

            if ((norm < tolerance .and. k > 1) .or. k > max_iterations) then
                exit
            end if

            ! Initialize the LHS matrix
            LHS = 0.0

            ! Populate LHS matrix
            do i = 2, N
                ! Momentum, Velocity
                LHS(i, i - 1) = LHS(i, i - 1) + 0.25*(-2.0*crossSectionLength(i - 1)*velocity(i - 1) - &
                            2.0*crossSectionLength(i)*velocity(i - 1) - &
                            crossSectionLength(i)*velocity(i) - crossSectionLength(i - 1)*velocity(i))
                LHS(i, i) = LHS(i, i) + crossSectionLength(i)*dx/tau + &
                            0.25*(crossSectionLength(i + 1)*velocity(i + 1) + &
                                    crossSectionLength(i)*velocity(i + 1) + &
                                    2.0*crossSectionLength(i + 1)*velocity(i) + &
                                    2.0*crossSectionLength(i)*velocity(i) - &
                                    crossSectionLength(i)*velocity(i - 1) - crossSectionLength(i - 1)*velocity(i - 1))
                LHS(i, i + 1) = LHS(i, i + 1) + 0.25*(crossSectionLength(i + 1)*velocity(i) + &
                            crossSectionLength(i)*velocity(i))

                ! Momentum, Pressure
                LHS(i, N + 1 + i - 1) = LHS(i, N + 1 + i - 1) - 0.25*crossSectionLength(i - 1) - &
                                                0.25*crossSectionLength(i)
                LHS(i, N + 1 + i) = LHS(i, N + 1 + i) + 0.25*crossSectionLength(i - 1) - &
                                            0.25*crossSectionLength(i + 1)
                LHS(i, N + 1 + i + 1) = LHS(i, N + 1 + i + 1) + 0.25*crossSectionLength(i) + &
                                                0.25*crossSectionLength(i + 1)
                ! Continuity, Velocity
                LHS(i + N + 1, i - 1) = LHS(i + N + 1, i - 1) - 0.25*crossSectionLength(i - 1) - &
                                0.25*crossSectionLength(i)
                LHS(i + N + 1, i) = LHS(i + N + 1, i) - 0.25*crossSectionLength(i - 1) + &
                                0.25*crossSectionLength(i + 1)
                LHS(i + N + 1, i + 1) = LHS(i + N + 1, i + 1) + 0.25*crossSectionLength(i) + &
                                0.25*crossSectionLength(i + 1)

                ! Continuity, Pressure
                LHS(i + N + 1, N + 1 + i - 1) = LHS(i + N + 1, N + 1 + i - 1) - alpha
                LHS(i + N + 1, N + 1 + i) = LHS(i + N + 1, N + 1 + i) + 2.0*alpha
                LHS(i + N + 1, N + 1 + i + 1) = LHS(i + N + 1, N + 1 + i + 1) - alpha
            end do

            ! Boundary conditions in LHS
            ! Velocity Inlet is prescribed
            LHS(1, 1) = 1.0
            ! Pressure Inlet is linearly interpolated
            LHS(N + 2, N + 2) = 1.0
            LHS(N + 2, N + 3) = -2.0
            LHS(N + 2, N + 4) = 1.0
            ! Velocity Outlet is linearly interpolated
            LHS(N + 1, N + 1) = 1.0
            LHS(N + 1, N) = -2.0
            LHS(N + 1, N - 1) = 1.0
            ! Pressure Outlet is Non-Reflecting
            LHS(2*N + 2, 2*N + 2) = 1.0
            LHS(2*N + 2, N + 1) = -(sqrt(c_mk2 - pressure_old(N + 1)/2.0) - (velocity(N + 1) - velocity_old(N + 1))/4.0)

            call dgesv(nlhs, nrhs, LHS, nlhs, ipiv, Res, nlhs, info)
            if (info /= 0) then
                write(*, *) "Linear Solver not converged!, Info: ", info
            end if

            ! Update velocity and pressure
            do i = 1, N + 1
                velocity(i) = velocity(i) + Res(i)
                pressure(i) = pressure(i) + Res(i + N + 1)
            end do
        end do

        ! Deallocate arrays
        deallocate(Res, LHS, ipiv)
        
    end subroutine fluid_compute_solution
end module FluidComputeSolution
