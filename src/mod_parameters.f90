module mod_parameters
  ! Reads para.in (Fortran NAMELIST) and holds all global
  ! simulation parameters. See para.in for the meaning of
  ! every entry.
  use mod_kinds
  implicit none

  integer(i4) :: N_subdivision
  real(dp)    :: R_apical, R_basal

  real(dp)    :: K_V, K_A, K_P, Lambda_line
  real(dp)    :: K_A_bas, K_P_bas
  real(dp)    :: V0_scale, A0_scale, A0_bas_scale

  real(dp)    :: zeta_friction, kBT, dt
  integer(i4) :: n_steps, it_dumps

  real(dp)    :: L_T1_threshold, A_T2_threshold
  integer(i4) :: it_topology_check
  logical     :: division_enable
  integer(i4) :: it_division_check
  real(dp)    :: V_division_threshold

  logical     :: T4_enable
  real(dp)    :: V_extrusion_threshold
  real(dp)    :: T4_direction_deadband

  logical     :: defect_enable
  real(dp)    :: defect_fraction
  real(dp)    :: Lambda_line_defect

  integer(i4) :: random_seed
  real(dp)    :: capacity_growth_factor

contains

  subroutine read_parameters(fname)
    character(len=*), intent(in) :: fname
    integer :: iun, ios
    namelist /SIMPARAMS/ N_subdivision, R_apical, R_basal, &
         K_V, K_A, K_P, Lambda_line, K_A_bas, K_P_bas, &
         V0_scale, A0_scale, A0_bas_scale, &
         zeta_friction, kBT, dt, n_steps, it_dumps, &
         L_T1_threshold, A_T2_threshold, it_topology_check, &
         division_enable, it_division_check, V_division_threshold, &
         T4_enable, V_extrusion_threshold, T4_direction_deadband, &
         defect_enable, defect_fraction, Lambda_line_defect, &
         random_seed, capacity_growth_factor

    ! default for the one optional/newer key, in case an older
    ! para.in without it is supplied
    capacity_growth_factor = 3.0_dp

    ! Basal-face moduli/target-scale default to a NEGATIVE sentinel
    ! here; if para.in does not set them explicitly, they are fixed up
    ! to equal the corresponding apical value right after the namelist
    ! read below (so an older para.in with no basal keys at all still
    ! runs, with the basal face behaving exactly like the apical one).
    K_A_bas      = -1.0_dp
    K_P_bas      = -1.0_dp
    A0_bas_scale = -1.0_dp

    ! T4_direction_deadband is the one optional T4 key (T4_enable and
    ! V_extrusion_threshold are required, like division_enable/
    ! V_division_threshold -- there is no sensible silent default for
    ! "should this feature be on" or "how compressed is crowded").
    T4_direction_deadband = -1.0_dp

    open(newunit=iun, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: cannot open parameter file: ', trim(fname)
      stop 1
    end if
    read(iun, nml=SIMPARAMS, iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: could not parse &SIMPARAMS namelist in ', trim(fname)
      stop 1
    end if
    close(iun)

    ! ---- basal-face default fix-up (see sentinel comment above) ----
    if (K_A_bas      < 0.0_dp) K_A_bas      = K_A
    if (K_P_bas      < 0.0_dp) K_P_bas      = K_P
    if (A0_bas_scale < 0.0_dp) A0_bas_scale = A0_scale
    if (T4_direction_deadband < 0.0_dp) T4_direction_deadband = 0.05_dp

    ! ---- sanity checks on the parameters themselves ----
    if (N_subdivision < 0 .or. N_subdivision > 6) then
      write(*,*) 'ERROR: N_subdivision out of sane range [0,6]: ', N_subdivision
      stop 1
    end if
    if (R_basal <= 0.0_dp .or. R_apical <= R_basal) then
      write(*,*) 'ERROR: need 0 < R_basal < R_apical. Got:', R_basal, R_apical
      stop 1
    end if
    if (dt <= 0.0_dp .or. n_steps <= 0) then
      write(*,*) 'ERROR: dt and n_steps must be positive.'
      stop 1
    end if
    if (zeta_friction <= 0.0_dp) then
      write(*,*) 'ERROR: zeta_friction must be positive.'
      stop 1
    end if
    if (capacity_growth_factor < 1.0_dp) then
      write(*,*) 'ERROR: capacity_growth_factor must be >= 1.0.'
      stop 1
    end if
    if (defect_enable .and. (defect_fraction <= 0.0_dp .or. defect_fraction > 1.0_dp)) then
      write(*,*) 'ERROR: defect_fraction must be in (0,1] when defect_enable is true. Got:', defect_fraction
      stop 1
    end if

    write(*,'(A)') '--- parameters read successfully ---'
    write(*,nml=SIMPARAMS)
    write(*,'(A)') '-------------------------------------'
  end subroutine read_parameters

end module mod_parameters
