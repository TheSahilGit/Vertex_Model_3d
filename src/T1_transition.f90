module mod_T1
  ! T1 ("neighbour exchange") junctional rearrangement.
  !
  ! Every interior vertex column is trivalent (exactly 3 cells meet
  ! there). Consider an apical edge (j,k) shared by cells A and B; let
  ! C be the third cell at j (other than A,B) and D the third cell at
  ! k. When the edge shrinks below L_T1_threshold we perform the
  ! standard T1 flip:
  !   - A loses vertex k (keeps j)         -> A shrinks by one side
  !   - B loses vertex j (keeps k)         -> B shrinks by one side
  !   - C gains k, inserted adjacent to j, on the side of j that used
  !     to border B                        -> C grows by one side
  !   - D gains j, inserted adjacent to k, on the side of k that used
  !     to border A                        -> D grows by one side
  ! so that after the flip the new edge (j,k) is shared by C and D
  ! instead of A and B, while both j and k remain trivalent points
  ! (j: {A,C,D}, k: {B,C,D}). The two vertices are then physically
  ! displaced apart along the direction perpendicular to the old edge
  ! (within the local tangent plane) so the new edge can re-open.
  !
  ! Same-side rule derivation: at the instant of the flip the 4 cells
  ! meet at a single 4-fold point with cyclic order (A,C,B,D); when it
  ! re-splits, the two new points separate into {A,C,D} and {B,C,D}.
  ! This is why C keeps its bond to A (gains k on the B-side) and D
  ! keeps its bond to B (gains j on the A-side).
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  use mod_topology
  implicit none

contains

  subroutine attempt_T1_transitions(n_done)
    integer(i4), intent(out) :: n_done
    logical, allocatable :: touched(:)
    integer(i4) :: ic, k, kp, n, j1, j2
    real(dp) :: elen

    allocate(touched(n_vert_cap))
    touched = .false.
    n_done = 0

    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      n = cells(ic)%n
      do k = 1, n
        kp = mod(k, n) + 1
        j1 = cells(ic)%vlist(k)
        j2 = cells(ic)%vlist(kp)
        if (touched(j1) .or. touched(j2)) cycle
        elen = norm3(r_api(:, j1) - r_api(:, j2))
        if (elen < L_T1_threshold) then
          call try_T1_edge(j1, j2, touched)
          if (touched(j1)) then
            n_done = n_done + 1
            ! cell ic is always one of the two edge-owning cells, so a
            ! successful flip just changed cells(ic)%n/%vlist under us
            ! (ic loses one vertex): the cached n/vlist for this cell
            ! are now stale, so stop scanning its (old) edge list here.
            exit
          end if
        end if
      end do
    end do

    deallocate(touched)
  end subroutine attempt_T1_transitions

  !------------------------------------------------------------------
  subroutine try_T1_edge(j, k, touched)
    integer(i4), intent(in)    :: j, k
    logical,     intent(inout) :: touched(:)

    integer(i4) :: A, B, C, D
    integer(i4) :: inc_j(8), inc_k(8), cnt_j, cnt_k, i
    integer(i4) :: afterC, afterD
    logical :: ok
    real(dp) :: mid_api(3), mid_bas(3), nrm(3), eold(3), tperp(3), reopen
    real(dp) :: nrm_len

    call edge_cells(j, k, A, B)
    if (A == 0 .or. B == 0) return

    call vertex_incidence(j, inc_j, cnt_j)
    call vertex_incidence(k, inc_k, cnt_k)
    if (cnt_j /= 3 .or. cnt_k /= 3) return   ! non-generic vertex; skip defensively

    C = 0
    do i = 1, cnt_j
      if (inc_j(i) /= A .and. inc_j(i) /= B) C = inc_j(i)
    end do
    D = 0
    do i = 1, cnt_k
      if (inc_k(i) /= A .and. inc_k(i) /= B) D = inc_k(i)
    end do
    if (C == 0 .or. D == 0 .or. C == D) return

    if (cells(A)%n <= 3 .or. cells(B)%n <= 3) return          ! would violate n>=3
    if (cells(C)%n + 1 > MAX_SIDES .or. cells(D)%n + 1 > MAX_SIDES) return

    call find_insert_side(C, j, B, afterC)
    if (afterC < 0) return
    call find_insert_side(D, k, A, afterD)
    if (afterD < 0) return

    ! ---- mutate topology (both insertions pre-checked to fit) ----
    call ring_insert_after(C, afterC, k, ok)
    if (.not. ok) return
    call ring_insert_after(D, afterD, j, ok)
    if (.not. ok) return
    call ring_remove(A, k, ok)
    call ring_remove(B, j, ok)

    ! ---- reposition j,k: separate perpendicular to the old edge ----
    eold = r_api(:, k) - r_api(:, j)
    mid_api = 0.5_dp * (r_api(:, j) + r_api(:, k))
    mid_bas = 0.5_dp * (r_bas(:, j) + r_bas(:, k))
    nrm = mid_api
    nrm_len = norm3(nrm)
    if (nrm_len > 1.0e-14_dp) nrm = nrm / nrm_len
    tperp = cross(nrm, eold)
    if (norm3(tperp) > 1.0e-14_dp) then
      tperp = tperp / norm3(tperp)
    else
      tperp = [1.0_dp, 0.0_dp, 0.0_dp]
    end if
    reopen = 1.5_dp * L_T1_threshold

    r_api(:, j) = mid_api - 0.5_dp * reopen * tperp
    r_api(:, k) = mid_api + 0.5_dp * reopen * tperp
    r_bas(:, j) = mid_bas - 0.5_dp * reopen * tperp
    r_bas(:, k) = mid_bas + 0.5_dp * reopen * tperp

    touched(j) = .true.
    touched(k) = .true.
  end subroutine try_T1_edge

  !------------------------------------------------------------------
  ! In icell's ring, vid has two neighbours (prev,next). Returns
  ! after_id such that ring_insert_after(icell, after_id, newv) places
  ! newv adjacent to vid on the side whose edge borders other_cell.
  ! Returns after_id = -1 if neither side matches (should not happen
  ! on a consistent mesh; caller aborts the T1 defensively).
  subroutine find_insert_side(icell, vid, other_cell, after_id)
    integer(i4), intent(in)  :: icell, vid, other_cell
    integer(i4), intent(out) :: after_id
    integer(i4) :: pos, n, prev_v, next_v

    after_id = -1
    pos = ring_find(icell, vid)
    if (pos == 0) return
    n = cells(icell)%n
    next_v = cells(icell)%vlist(mod(pos, n) + 1)
    prev_v = cells(icell)%vlist(mod(pos - 2 + n, n) + 1)

    if (ring_has_edge(other_cell, vid, next_v)) then
      after_id = vid
    else if (ring_has_edge(other_cell, vid, prev_v)) then
      after_id = prev_v
    end if
  end subroutine find_insert_side

end module mod_T1
