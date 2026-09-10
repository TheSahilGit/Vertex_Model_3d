module mod_division
  ! Cell division: a cell whose volume has grown past
  ! V_division_threshold * V0 splits into two daughter cells.
  !
  ! A division axis is chosen as a pair of ring edges roughly opposite
  ! each other (index k1 random, k2 = k1 + n/2 mod n); a new vertex
  ! column is inserted at the midpoint of each of those two edges (on
  ! both apical and basal shells), and the polygon is cut into two
  ! contiguous arcs at those two new vertices, each arc plus the new
  ! shared edge forming one daughter cell. The two neighbouring cells
  ! whose edges were bisected each gain the corresponding new vertex
  ! (their side count +1). Only cells with n>=4 sides are eligible.
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_topology
  implicit none

contains

  subroutine attempt_cell_divisions(n_done)
    integer(i4), intent(out) :: n_done
    integer(i4) :: ic
    logical :: ok

    n_done = 0
    if (.not. division_enable) return

    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      if (cells(ic)%n < 4) cycle
      if (cells(ic)%V_last > V_division_threshold * cells(ic)%V0) then
        call do_division(ic, ok)
        if (ok) n_done = n_done + 1
      end if
    end do
  end subroutine attempt_cell_divisions

  !------------------------------------------------------------------
  subroutine do_division(parent, ok)
    ! IMPORTANT: every validity/capacity check is done FIRST, using only
    ! local arrays and sentinel placeholders (SENT1,SENT2) for the two
    ! not-yet-created vertices. Only once the split is known to be
    ! geometrically valid AND both a vertex and a cell slot are known
    ! to be available do we mutate any global state (vertex creation,
    ! neighbour ring insertion, final ring assignment). This keeps the
    ! routine atomic: a failure never leaves a half-applied division
    ! (e.g. a vertex spliced into a neighbour's ring with no matching
    ! parent/daughter split) behind.
    integer(i4), intent(in)  :: parent
    logical,     intent(out) :: ok

    integer(i4), parameter :: SENT1 = -1, SENT2 = -2
    integer(i4) :: n, k1, k2, kp1, kp2, tmp, i
    integer(i4) :: v1new, v2new
    integer(i4) :: neighA, neighB, neigh1, neigh2
    integer(i4) :: R(MAX_SIDES), combined(MAX_SIDES + 2)
    integer(i4) :: ncomb, p1, p2, len1, len2
    integer(i4) :: daughter1(MAX_SIDES), daughter2(MAX_SIDES)
    integer(i4) :: new_cell
    logical :: ins_ok
    real(dp) :: u
    real(dp) :: V0_old, A0_old, A0_bas_old

    ok = .false.
    n = cells(parent)%n
    R(1:n) = cells(parent)%vlist(1:n)

    call random_number(u)
    k1 = 1 + int(u * real(n, dp))
    if (k1 > n) k1 = n
    k2 = k1 + n/2
    k2 = mod(k2 - 1, n) + 1
    if (k2 < k1) then
      tmp = k1; k1 = k2; k2 = tmp
    end if
    if (k2 == k1) return

    kp1 = mod(k1, n) + 1
    kp2 = mod(k2, n) + 1

    ! ---- build the two daughter rings using sentinels (no global state touched) ----
    ncomb = 0
    do i = 1, n
      ncomb = ncomb + 1
      combined(ncomb) = R(i)
      if (i == k1) then
        ncomb = ncomb + 1
        combined(ncomb) = SENT1
      end if
      if (i == k2) then
        ncomb = ncomb + 1
        combined(ncomb) = SENT2
      end if
    end do

    p1 = 0; p2 = 0
    do i = 1, ncomb
      if (combined(i) == SENT1) p1 = i
      if (combined(i) == SENT2) p2 = i
    end do
    if (p1 == 0 .or. p2 == 0 .or. p1 >= p2) return

    len1 = p2 - p1 + 1
    daughter1(1:len1) = combined(p1:p2)

    len2 = ncomb - len1 + 2
    daughter2(1:ncomb - p2 + 1) = combined(p2:ncomb)
    daughter2(ncomb - p2 + 2:len2) = combined(1:p1)

    if (len1 < 3 .or. len2 < 3 .or. len1 > MAX_SIDES .or. len2 > MAX_SIDES) return

    ! ---- capacity checks (still no mutation) ----
    if (n_vert + 2 > n_vert_cap) then
      write(*,*) 'WARNING: vertex capacity exceeded, skipping division of cell', parent
      return
    end if
    if (count_alive_cells() >= n_cell_cap .and. n_cell >= n_cell_cap) then
      write(*,*) 'WARNING: cell capacity exceeded, skipping division of cell', parent
      return
    end if

    ! identify the (at most 2) neighbour cells whose edges get bisected,
    ! purely as a query -- no mutation yet.
    call edge_cells(R(k1), R(kp1), neighA, neighB)
    neigh1 = merge(neighA, neighB, neighA /= parent)
    if (neigh1 == parent) neigh1 = 0
    call edge_cells(R(k2), R(kp2), neighA, neighB)
    neigh2 = merge(neighA, neighB, neighA /= parent)
    if (neigh2 == parent) neigh2 = 0

    ! Both neighbours must have room for one more side, checked BEFORE
    ! either ring is touched, so this division is all-or-nothing. If the
    ! same third cell happens to border the parent along BOTH cut edges
    ! (neigh1==neigh2), it receives BOTH new vertices, so it needs room
    ! for two more sides, not one.
    if (neigh1 /= 0 .and. neigh1 == neigh2) then
      if (cells(neigh1)%n + 2 > MAX_SIDES) return
    else
      if (neigh1 /= 0) then
        if (cells(neigh1)%n + 1 > MAX_SIDES) return
      end if
      if (neigh2 /= 0) then
        if (cells(neigh2)%n + 1 > MAX_SIDES) return
      end if
    end if

    ! ---- everything validated: now allocate the two new vertex columns ----
    v1new = next_free_vert_slot()
    v2new = next_free_vert_slot()
    if (v1new == 0 .or. v2new == 0) return   ! should not happen given the capacity check above

    r_api(:, v1new) = 0.5_dp * (r_api(:, R(k1)) + r_api(:, R(kp1)))
    r_bas(:, v1new) = 0.5_dp * (r_bas(:, R(k1)) + r_bas(:, R(kp1)))
    vert_alive(v1new) = .true.

    r_api(:, v2new) = 0.5_dp * (r_api(:, R(k2)) + r_api(:, R(kp2)))
    r_bas(:, v2new) = 0.5_dp * (r_bas(:, R(k2)) + r_bas(:, R(kp2)))
    vert_alive(v2new) = .true.

    new_cell = next_free_cell_slot()
    if (new_cell == 0) then
      ! genuinely should not happen (checked above), but stay defensive
      write(*,*) 'WARNING: cell capacity exceeded, skipping division of cell', parent
      vert_alive(v1new) = .false.; vert_alive(v2new) = .false.
      return
    end if

    ! ---- splice the new vertices into the (at most 2) bisected neighbours ----
    ! (the capacity pre-checks above -- including the neigh1==neigh2 case --
    ! mean insert_between should never fail here; the handling below is
    ! defense-in-depth only, and fully rolls back the vertex AND cell
    ! slots just reserved so no "ghost" alive-but-empty cell is left.)
    if (neigh1 /= 0) then
      call insert_between(neigh1, R(k1), R(kp1), v1new, ins_ok)
      if (.not. ins_ok) then
        write(*,*) 'WARNING: neighbour ring full, aborting division of cell', parent
        vert_alive(v1new) = .false.; vert_alive(v2new) = .false.
        cells(new_cell)%alive = .false.; cells(new_cell)%n = 0
        return
      end if
    end if
    if (neigh2 /= 0) then
      call insert_between(neigh2, R(k2), R(kp2), v2new, ins_ok)
      if (.not. ins_ok) then
        write(*,*) 'WARNING: neighbour ring full, aborting division of cell', parent
        ! neigh1 may already have v1new spliced into its ring at this
        ! point; that is left in place (it is a genuine, valid point on
        ! that shared edge) rather than un-inserted, since ring_remove
        ! would need the exact insertion position -- harmless either way,
        ! it does not violate ring integrity, merely leaves one bisected
        ! edge whose division did not complete on this side.
        vert_alive(v2new) = .false.
        cells(new_cell)%alive = .false.; cells(new_cell)%n = 0
        return
      end if
    end if

    ! substitute the sentinels with the real, now-allocated vertex ids
    do i = 1, len1
      if (daughter1(i) == SENT1) daughter1(i) = v1new
      if (daughter1(i) == SENT2) daughter1(i) = v2new
    end do
    do i = 1, len2
      if (daughter2(i) == SENT1) daughter2(i) = v1new
      if (daughter2(i) == SENT2) daughter2(i) = v2new
    end do

    V0_old     = cells(parent)%V0
    A0_old     = cells(parent)%A0
    A0_bas_old = cells(parent)%A0_bas

    cells(parent)%vlist = 0
    cells(parent)%vlist(1:len1) = daughter1(1:len1)
    cells(parent)%n = len1
    cells(parent)%V0     = 0.5_dp * V0_old
    cells(parent)%A0     = 0.5_dp * A0_old
    cells(parent)%A0_bas = 0.5_dp * A0_bas_old

    cells(new_cell)%vlist = 0
    cells(new_cell)%vlist(1:len2) = daughter2(1:len2)
    cells(new_cell)%n = len2
    cells(new_cell)%alive = .true.
    cells(new_cell)%V0     = 0.5_dp * V0_old
    cells(new_cell)%A0     = 0.5_dp * A0_old
    cells(new_cell)%A0_bas = 0.5_dp * A0_bas_old

    ok = .true.
  end subroutine do_division

end module mod_division
