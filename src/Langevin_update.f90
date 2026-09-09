module mod_langevin
  ! Overdamped Langevin integrator:
  !     zeta * dr/dt = -dE/dr + xi(t),   <xi_a(t) xi_b(t')> = 2*zeta*kBT*delta_ab*delta(t-t')
  ! discretised (Euler-Maruyama) as
  !     r(t+dt) = r(t) + (dt/zeta) * F(t) + sqrt(2*kBT*dt/zeta) * N(0,1)
  ! applied independently to every apical and basal vertex position.
  use mod_kinds
  use mod_parameters
  use mod_data
  implicit none

contains

  subroutine init_rng(seed)
    integer(i4), intent(in) :: seed
    integer :: n, i
    integer, allocatable :: seed_arr(:)
    call random_seed(size=n)
    allocate(seed_arr(n))
    do i = 1, n
      seed_arr(i) = seed + i * 104729
    end do
    call random_seed(put=seed_arr)
    deallocate(seed_arr)
  end subroutine init_rng

  real(dp) function gaussian_rand() result(g)
    real(dp) :: u1, u2
    real(dp), parameter :: two_pi = 6.283185307179586_dp
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, 1.0e-300_dp)
    g = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
  end function gaussian_rand

  subroutine langevin_step()
    real(dp) :: coef, noise_amp
    integer(i4) :: i, d

    coef = dt / zeta_friction
    noise_amp = sqrt(2.0_dp * kBT * dt / zeta_friction)

    do i = 1, n_vert
      if (.not. vert_alive(i)) cycle
      do d = 1, 3
        r_api(d, i) = r_api(d, i) + coef * f_api(d, i)
        r_bas(d, i) = r_bas(d, i) + coef * f_bas(d, i)
        if (noise_amp > 0.0_dp) then
          r_api(d, i) = r_api(d, i) + noise_amp * gaussian_rand()
          r_bas(d, i) = r_bas(d, i) + noise_amp * gaussian_rand()
        end if
      end do
    end do
  end subroutine langevin_step

end module mod_langevin
