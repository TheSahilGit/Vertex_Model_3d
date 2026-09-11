module mod_topology
  ! Low-level combinatorial helpers on the cell-ring data structure
  ! (mod_data % cells(:) % vlist). Used by T1_transition, T2_transition
  ! and Cell_Division to perform the actual "ring surgery". Kept
  ! independent of geometry (mod_geometry) on purpose.
  use mod_kinds
  use mod_data
  implicit none

contains

  !------------------------------------------------------------------
  integer(i4) function ring_find(icell, vid) result(pos)
    integer(i4), intent(in) :: icell, vid
    integer(i4) :: k
    pos = 0
    do k = 1, cells(icell)%n
      if (cells(icell)%vlist(k) == vid) then
        pos = k
        return
      end if
    end do
  end function ring_find

  !------------------------------------------------------------------
  logical function ring_has_edge(icell, v1, v2) result(has)
    integer(i4), intent(in) :: icell, v1, v2
    integer(i4) :: k, n, a, b
    has = .false.
    n = cells(icell)%n
    do k = 1, n
      a = cells(icell)%vlist(k)
      b = cells(icell)%vlist(mod(k, n) + 1)
      if ((a == v1 .and. b == v2) .or. (a == v2 .and. b == v1)) then
        has = .true.
        return
      end if
    end do
  end function ring_has_edge

  !------------------------------------------------------------------
  ! Find the (up to two) alive cells whose ring contains the ordered
  ! (cyclic, either direction) adjacent pair (v1,v2). Returns 0 for a
  ! slot that is not found (should not normally happen on a closed
  ! manifold mesh -- callers should treat 0 as "abort this event").
  subroutine edge_cells(v1, v2, cA, cB)
    integer(i4), intent(in)  :: v1, v2
    integer(i4), intent(out) :: cA, cB
    integer(i4) :: i, nfound
    cA = 0; cB = 0; nfound = 0
    do i = 1, n_cell
      if (.not. cells(i)%alive) cycle
      if (ring_has_edge(i, v1, v2)) then
        nfound = nfound + 1
        if (nfound == 1) then
          cA = i
        else if (nfound == 2) then
          cB = i
        end if
      end if
    end do
  end subroutine edge_cells

  !------------------------------------------------------------------
  ! All alive cells containing vertex column vid (generically 3 of them).
  subroutine vertex_incidence(vid, list, count)
    integer(i4), intent(in)  :: vid
    integer(i4), intent(out) :: list(:)
    integer(i4), intent(out) :: count
    integer(i4) :: i
    count = 0
    do i = 1, n_cell
      if (.not. cells(i)%alive) cycle
      if (ring_find(i, vid) > 0) then
        count = count + 1
        if (count <= size(list)) list(count) = i
      end if
    end do
  end subroutine vertex_incidence

  !------------------------------------------------------------------
  ! Cheap O(n_cell)-TOTAL (not per-edge) "which alive cells touch each
  ! vertex" cache, built once and reused for many neighbour_across_edge
  ! lookups -- used by Force.f90's defect_boundary_tension feature.
  ! Calling edge_cells (O(n_cell) EACH) once per lateral face would make
  ! compute_forces (called every step) effectively O(n_cell^2); building
  ! this cache once per compute_forces call instead keeps it the same
  ! O(n_cell*MAX_SIDES) order as the rest of that routine.
  ! vic(1:3, v) holds up to 3 owning cell ids (0-padded); a vertex
  ! touched by more than 3 alive cells (should not happen on this
  ! trivalent mesh) simply has the extras silently dropped.
  subroutine build_vertex_incidence(vic, vic_count)
    integer(i4), intent(out) :: vic(:,:)
    integer(i4), intent(out) :: vic_count(:)
    integer(i4) :: ic, k, v

    vic = 0
    vic_count = 0
    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      do k = 1, cells(ic)%n
        v = cells(ic)%vlist(k)
        vic_count(v) = vic_count(v) + 1
        if (vic_count(v) <= size(vic, 1)) vic(vic_count(v), v) = ic
      end do
    end do
  end subroutine build_vertex_incidence

  !------------------------------------------------------------------
  ! Given icell's own edge (v1,v2) and the cache above, return the
  ! OTHER cell owning that edge (0 if not found -- should not happen on
  ! a closed manifold mesh).
  integer(i4) function neighbor_across_edge(icell, v1, v2, vic, vic_count) result(jc)
    integer(i4), intent(in) :: icell, v1, v2
    integer(i4), intent(in) :: vic(:,:), vic_count(:)
    integer(i4) :: a, b, cand

    jc = 0
    do a = 1, min(vic_count(v1), size(vic, 1))
      cand = vic(a, v1)
      if (cand == icell .or. cand == 0) cycle
      do b = 1, min(vic_count(v2), size(vic, 1))
        if (vic(b, v2) == cand) then
          jc = cand
          return
        end if
      end do
    end do
  end function neighbor_across_edge

  !------------------------------------------------------------------
  subroutine ring_remove(icell, vid, ok)
    integer(i4), intent(in)  :: icell, vid
    logical,     intent(out) :: ok
    integer(i4) :: pos, k, n
    pos = ring_find(icell, vid)
    if (pos == 0) then
      ok = .false.
      return
    end if
    n = cells(icell)%n
    do k = pos, n - 1
      cells(icell)%vlist(k) = cells(icell)%vlist(k + 1)
    end do
    cells(icell)%vlist(n) = 0
    cells(icell)%n = n - 1
    ok = (cells(icell)%n >= 3)
  end subroutine ring_remove

  !------------------------------------------------------------------
  subroutine ring_insert_after(icell, after_vid, new_vid, ok)
    integer(i4), intent(in)  :: icell, after_vid, new_vid
    logical,     intent(out) :: ok
    integer(i4) :: pos, k, n
    pos = ring_find(icell, after_vid)
    if (pos == 0) then
      ok = .false.
      return
    end if
    n = cells(icell)%n
    if (n + 1 > MAX_SIDES) then
      ok = .false.
      return
    end if
    do k = n, pos + 1, -1
      cells(icell)%vlist(k + 1) = cells(icell)%vlist(k)
    end do
    cells(icell)%vlist(pos + 1) = new_vid
    cells(icell)%n = n + 1
    ok = .true.
  end subroutine ring_insert_after

  !------------------------------------------------------------------
  ! Replace every occurrence of old_vid in the ring by new_vid (size
  ! unchanged). Caller should follow with collapse_duplicates.
  subroutine ring_replace(icell, old_vid, new_vid)
    integer(i4), intent(in) :: icell, old_vid, new_vid
    integer(i4) :: k
    do k = 1, cells(icell)%n
      if (cells(icell)%vlist(k) == old_vid) cells(icell)%vlist(k) = new_vid
    end do
  end subroutine ring_replace

  !------------------------------------------------------------------
  ! Collapse consecutive (cyclically) duplicate ring entries, e.g. after
  ! a T2 vertex merge maps two ring-adjacent ids onto the same value.
  subroutine collapse_duplicates(icell)
    integer(i4), intent(in) :: icell
    integer(i4) :: n, k, m, tmp(MAX_SIDES)
    n = cells(icell)%n
    if (n <= 1) return
    tmp(1) = cells(icell)%vlist(1)
    m = 1
    do k = 2, n
      if (cells(icell)%vlist(k) /= tmp(m)) then
        m = m + 1
        tmp(m) = cells(icell)%vlist(k)
      end if
    end do
    if (m >= 2) then
      if (tmp(m) == tmp(1)) m = m - 1
    end if
    cells(icell)%vlist = 0
    cells(icell)%vlist(1:m) = tmp(1:m)
    cells(icell)%n = m
  end subroutine collapse_duplicates

  !------------------------------------------------------------------
  ! Insert newv into icell's ring so that it sits between the two
  ! (already ring-adjacent) vertices va,vb, regardless of which order
  ! they currently appear in (two cells sharing an edge walk it in
  ! opposite directions in their own rings).
  subroutine insert_between(icell, va, vb, newv, ok)
    integer(i4), intent(in)  :: icell, va, vb, newv
    logical,     intent(out) :: ok
    integer(i4) :: posa, n, nxt

    ok = .false.
    posa = ring_find(icell, va)
    if (posa == 0) return
    n = cells(icell)%n
    nxt = cells(icell)%vlist(mod(posa, n) + 1)
    if (nxt == vb) then
      call ring_insert_after(icell, va, newv, ok)
    else
      call ring_insert_after(icell, vb, newv, ok)
    end if
  end subroutine insert_between

  !------------------------------------------------------------------
  ! Both slot-allocator functions mark the returned slot as alive/used
  ! BEFORE returning (not left to the caller): they may be called
  ! several times back-to-back (e.g. twice per cell division, for the
  ! two new vertex columns) and must never hand out the same free slot
  ! twice. The caller still owns setting the slot's real content
  ! (position, ring, V0/A0, ...); "alive" merely reserves it.
  integer(i4) function next_free_cell_slot() result(id)
    integer(i4) :: i
    id = 0
    do i = 1, n_cell
      if (.not. cells(i)%alive) then
        id = i
        cells(id)%alive = .true.
        return
      end if
    end do
    if (n_cell < n_cell_cap) then
      n_cell = n_cell + 1
      id = n_cell
      cells(id)%alive = .true.
    end if
  end function next_free_cell_slot

  !------------------------------------------------------------------
  integer(i4) function next_free_vert_slot() result(id)
    integer(i4) :: i
    id = 0
    do i = 1, n_vert
      if (.not. vert_alive(i)) then
        id = i
        vert_alive(id) = .true.
        return
      end if
    end do
    if (n_vert < n_vert_cap) then
      n_vert = n_vert + 1
      id = n_vert
      vert_alive(id) = .true.
    end if
  end function next_free_vert_slot

  !------------------------------------------------------------------
  ! Euler-characteristic sanity check for the current alive mesh:
  ! V - E + F = 2 for a closed genus-0 surface. F = number of alive
  ! cells, V = number of alive vertex columns, E = number of alive
  ! edges (each counted once even though shared by 2 cells: we count
  ! Sum(n_sides)/2).
  subroutine euler_check(ok, V, E, F, chi)
    logical, intent(out)     :: ok
    integer(i4), intent(out) :: V, E, F
    real(dp), intent(out)    :: chi
    integer(i4) :: i, sum_sides

    V = count_alive_verts()
    F = 0
    sum_sides = 0
    do i = 1, n_cell
      if (cells(i)%alive) then
        F = F + 1
        sum_sides = sum_sides + cells(i)%n
      end if
    end do
    if (mod(sum_sides, 2) /= 0) then
      E = -1
      chi = -999.0_dp
      ok = .false.
      return
    end if
    E = sum_sides / 2
    chi = real(V - E + F, dp)
    ok = (abs(chi - 2.0_dp) < 1.0e-8_dp)
  end subroutine euler_check

end module mod_topology
