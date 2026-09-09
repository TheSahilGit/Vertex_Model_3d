module mod_sanity
  ! Automated self-tests run once at start-up (before/independently of
  ! the real mesh), plus a few runtime diagnostic checks called during
  ! the simulation. These exist because a wrong sign/winding convention
  ! in the geometry code would silently produce a plausible-looking but
  ! physically wrong simulation -- catching it here, against known
  ! analytic answers, is much cheaper than discovering it later.
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  use mod_topology
  implicit none

  real(dp), parameter :: pi_ = 3.14159265358979323846_dp

contains

  !=====================================================================
  ! Finite-difference check of tri_area_grad and tetra_vol_grad on a
  ! fixed, generic (non-degenerate, non-symmetric) triangle.
  subroutine test_geometry_gradients(ok)
    logical, intent(out) :: ok
    real(dp) :: p(3), q(3), r(3)
    real(dp) :: area, gp(3), gq(3), gr(3)
    real(dp) :: vol, hp(3), hq(3), hr(3)
    real(dp) :: fd(3), err
    real(dp), parameter :: h = 1.0e-6_dp, tol = 1.0e-6_dp

    ok = .true.
    p = [0.10_dp,  0.20_dp,  0.30_dp]
    q = [0.90_dp, -0.40_dp,  0.15_dp]
    r = [-0.30_dp, 0.65_dp, -0.50_dp]

    call tri_area_grad(p, q, r, area, gp, gq, gr)
    call fd_grad_area(p, q, r, 1, h, fd); err = norm3(fd - gp)
    if (err > tol) then; ok=.false.; write(*,*) 'FAIL tri_area_grad wrt p, err=', err; end if
    call fd_grad_area(p, q, r, 2, h, fd); err = norm3(fd - gq)
    if (err > tol) then; ok=.false.; write(*,*) 'FAIL tri_area_grad wrt q, err=', err; end if
    call fd_grad_area(p, q, r, 3, h, fd); err = norm3(fd - gr)
    if (err > tol) then; ok=.false.; write(*,*) 'FAIL tri_area_grad wrt r, err=', err; end if

    call tetra_vol_grad(p, q, r, vol, hp, hq, hr)
    call fd_grad_vol(p, q, r, 1, h, fd); err = norm3(fd - hp)
    if (err > tol) then; ok=.false.; write(*,*) 'FAIL tetra_vol_grad wrt p, err=', err; end if
    call fd_grad_vol(p, q, r, 2, h, fd); err = norm3(fd - hq)
    if (err > tol) then; ok=.false.; write(*,*) 'FAIL tetra_vol_grad wrt q, err=', err; end if
    call fd_grad_vol(p, q, r, 3, h, fd); err = norm3(fd - hr)
    if (err > tol) then; ok=.false.; write(*,*) 'FAIL tetra_vol_grad wrt r, err=', err; end if

    if (ok) write(*,'(A)') 'sanity: analytic area/volume gradients match finite differences. OK'

  contains
    subroutine fd_grad_area(p, q, r, which, h, fd)
      real(dp), intent(in) :: p(3), q(3), r(3), h
      integer, intent(in) :: which
      real(dp), intent(out) :: fd(3)
      real(dp) :: pp(3), pm(3), qq(3), qm(3), rr(3), rm(3)
      real(dp) :: aplus, aminus, dumv(3), dumv2(3), dumv3(3)
      integer :: d
      do d = 1, 3
        pp = p; qq = q; rr = r
        pm = p; qm = q; rm = r
        select case (which)
        case (1); pp(d) = pp(d) + h; pm(d) = pm(d) - h
        case (2); qq(d) = qq(d) + h; qm(d) = qm(d) - h
        case (3); rr(d) = rr(d) + h; rm(d) = rm(d) - h
        end select
        call tri_area_grad(pp, qq, rr, aplus, dumv, dumv2, dumv3)
        call tri_area_grad(pm, qm, rm, aminus, dumv, dumv2, dumv3)
        fd(d) = (aplus - aminus) / (2.0_dp * h)
      end do
    end subroutine fd_grad_area

    subroutine fd_grad_vol(p, q, r, which, h, fd)
      real(dp), intent(in) :: p(3), q(3), r(3), h
      integer, intent(in) :: which
      real(dp), intent(out) :: fd(3)
      real(dp) :: pp(3), pm(3), qq(3), qm(3), rr(3), rm(3)
      real(dp) :: vplus, vminus, dumv(3), dumv2(3), dumv3(3)
      integer :: d
      do d = 1, 3
        pp = p; qq = q; rr = r
        pm = p; qm = q; rm = r
        select case (which)
        case (1); pp(d) = pp(d) + h; pm(d) = pm(d) - h
        case (2); qq(d) = qq(d) + h; qm(d) = qm(d) - h
        case (3); rr(d) = rr(d) + h; rm(d) = rm(d) - h
        end select
        call tetra_vol_grad(pp, qq, rr, vplus, dumv, dumv2, dumv3)
        call tetra_vol_grad(pm, qm, rm, vminus, dumv, dumv2, dumv3)
        fd(d) = (vplus - vminus) / (2.0_dp * h)
      end do
    end subroutine fd_grad_vol
  end subroutine test_geometry_gradients

  !=====================================================================
  ! Build a single unit-cube "cell" (4 vertex columns: apical = top
  ! square z=1, basal = bottom square z=0) and check that the SAME
  ! production routine used by the real simulation (compute_forces)
  ! recovers volume=1 and apical area=1. This directly validates the
  ! face-winding / triangulation convention used throughout Force.f90.
  subroutine test_cube_volume(ok)
    use mod_force, only: compute_forces
    logical, intent(out) :: ok
    real(dp) :: energy, vol_total, area_total

    call allocate_data(1_i4, 4_i4, growth_factor=1.0_dp)
    n_cell = 1; n_vert = 4

    r_api(:,1) = [0.0_dp, 0.0_dp, 1.0_dp]
    r_api(:,2) = [1.0_dp, 0.0_dp, 1.0_dp]
    r_api(:,3) = [1.0_dp, 1.0_dp, 1.0_dp]
    r_api(:,4) = [0.0_dp, 1.0_dp, 1.0_dp]
    r_bas(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    r_bas(:,2) = [1.0_dp, 0.0_dp, 0.0_dp]
    r_bas(:,3) = [1.0_dp, 1.0_dp, 0.0_dp]
    r_bas(:,4) = [0.0_dp, 1.0_dp, 0.0_dp]
    vert_alive(1:4) = .true.

    cells(1)%n = 4
    cells(1)%vlist = 0
    cells(1)%vlist(1:4) = [1,2,3,4]
    cells(1)%alive = .true.
    cells(1)%V0 = 0.0_dp
    cells(1)%A0 = 0.0_dp

    call compute_forces(energy, vol_total, area_total)

    ok = (abs(cells(1)%V_last - 1.0_dp) < 1.0e-10_dp) .and. &
         (abs(cells(1)%A_last - 1.0_dp) < 1.0e-10_dp)

    if (ok) then
      write(*,'(A)') 'sanity: unit-cube volume/area test passed (V=1, Aapi=1). OK'
    else
      write(*,'(A,F12.8,A,F12.8)') 'FAIL unit-cube test: V_last=', cells(1)%V_last, &
           '  A_last=', cells(1)%A_last
    end if
  end subroutine test_cube_volume

  !=====================================================================
  subroutine global_shell_checks(vol_total, area_total)
    real(dp), intent(in) :: vol_total, area_total
    real(dp) :: expected_vol, expected_area, relerr_v, relerr_a

    expected_area = 4.0_dp * pi_ * R_apical**2
    expected_vol  = (4.0_dp/3.0_dp) * pi_ * (R_apical**3 - R_basal**3)
    relerr_a = abs(area_total - expected_area) / expected_area
    relerr_v = abs(vol_total  - expected_vol ) / expected_vol

    write(*,'(A,ES14.6,A,ES14.6,A,F6.3)') 'sanity: total apical area = ', area_total, &
         '  (ideal sphere = ', expected_area, ')  rel.err=', relerr_a
    write(*,'(A,ES14.6,A,ES14.6,A,F6.3)') 'sanity: total shell volume = ', vol_total, &
         '  (ideal shell  = ', expected_vol,  ')  rel.err=', relerr_v

    if (relerr_a > 0.5_dp .or. relerr_v > 0.5_dp) then
      write(*,'(A)') 'WARNING: total area/volume deviates >50% from the ideal shell -- check for a bug.'
    end if
  end subroutine global_shell_checks

  !=====================================================================
  ! Checks that every alive cell's ring vlist(1:n) contains n distinct,
  ! valid (1<=id<=n_vert, alive) vertex-column ids -- i.e. no stale
  ! zeros or dangling references. Prints the first violation found.
  subroutine ring_integrity_check(tag, ok)
    character(len=*), intent(in) :: tag
    logical, intent(out) :: ok
    integer(i4) :: ic, k, kk, n, vid
    ok = .true.
    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      n = cells(ic)%n
      if (n < 3 .or. n > MAX_SIDES) then
        write(*,'(A,A,A,I0,A,I0)') 'RING-CHECK[', tag, '] FAIL: cell ', ic, ' has bad n=', n
        ok = .false.; return
      end if
      do k = 1, n
        vid = cells(ic)%vlist(k)
        if (vid < 1 .or. vid > n_vert) then
          write(*,'(A,A,A,I0,A,I0,A,I0)') 'RING-CHECK[', tag, '] FAIL: cell ', ic, ' slot ', k, ' has id=', vid
          ok = .false.; return
        end if
        if (.not. vert_alive(vid)) then
          write(*,'(A,A,A,I0,A,I0,A,I0)') 'RING-CHECK[', tag, '] FAIL: cell ', ic, ' slot ', k, ' refers to dead vertex ', vid
          ok = .false.; return
        end if
        do kk = k+1, n
          if (cells(ic)%vlist(kk) == vid) then
            write(*,'(A,A,A,I0,A,I0)') 'RING-CHECK[', tag, '] FAIL: cell ', ic, ' has duplicate id ', vid
            ok = .false.; return
          end if
        end do
      end do
    end do
    if (ok) write(*,'(A,A,A)') 'sanity [', tag, ']: ring integrity OK.'
  end subroutine ring_integrity_check

  !=====================================================================
  subroutine topology_check(tag)
    character(len=*), intent(in) :: tag
    logical :: ok
    integer(i4) :: V, E, F
    real(dp) :: chi
    call euler_check(ok, V, E, F, chi)
    if (.not. ok) then
      write(*,'(A,A,A,I0,A,I0,A,I0,A,F8.3)') 'WARNING [', tag, ']: Euler characteristic check FAILED. V=', &
           V, ' E=', E, ' F=', F, ' chi=', chi
    else
      write(*,'(A,A,A,I0,A,I0,A,I0)') 'sanity [', tag, ']: Euler V-E+F=2 OK.  V=', V, ' E=', E, ' F=', F
    end if
  end subroutine topology_check

end module mod_sanity
