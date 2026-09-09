module mod_T2
  ! T2 transition: elimination of a triangular cell whose apical area
  ! has shrunk below A_T2_threshold (the standard vertex-model
  ! counterpart of cell extrusion/death). Only n==3 cells are eligible
  ! -- a cell reaches 3 sides through a sequence of T1 events before it
  ! can collapse, which is the usual convention.
  !
  ! The triangle's 3 vertex columns are merged into a single point
  ! (their mean, on both apical and basal shells); one of the three
  ! ids is kept as the surviving vertex and the other two are replaced
  ! by it EVERYWHERE (a global substitution over every other cell's
  ! ring), after which any resulting consecutive duplicate ring entries
  ! are collapsed. This shrinks each of the (up to three) neighbouring
  ! cells by exactly one side, which is the standard effect of a T2
  ! event -- EXCEPT when a neighbour itself only had 3 or 4 sides, in
  ! which case this can push it below 3 sides (a cascading collapse:
  ! e.g. two T2 "victim" triangles sharing an edge). Rather than leave
  ! such a degenerate cell behind, it is queued and forcibly merged
  ! down the same way (its remaining 2 or 1 vertices collapsed to a
  ! point and the cell removed), which may itself cascade further; the
  ! queue is drained until stable (bounded, so a pathological mesh
  ! cannot loop forever).
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  use mod_topology
  implicit none

contains

  subroutine attempt_T2_transitions(n_done)
    integer(i4), intent(out) :: n_done
    integer(i4) :: ic
    real(dp) :: area

    n_done = 0
    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      if (cells(ic)%n /= 3) cycle
      area = triangle_cell_area(ic)
      if (area < A_T2_threshold) then
        call do_T2(ic)
        n_done = n_done + 1
      end if
    end do
  end subroutine attempt_T2_transitions

  real(dp) function triangle_cell_area(icell) result(area)
    integer(i4), intent(in) :: icell
    integer(i4) :: j1, j2, j3
    real(dp) :: n(3)
    j1 = cells(icell)%vlist(1); j2 = cells(icell)%vlist(2); j3 = cells(icell)%vlist(3)
    n = cross(r_api(:, j2) - r_api(:, j1), r_api(:, j3) - r_api(:, j1))
    area = 0.5_dp * norm3(n)
  end function triangle_cell_area

  !------------------------------------------------------------------
  subroutine do_T2(icell)
    integer(i4), intent(in) :: icell
    integer(i4), parameter :: QCAP = 256
    integer(i4) :: queue(QCAP), qhead, qtail
    integer(i4) :: cur

    qhead = 1; qtail = 0

    qtail = qtail + 1; queue(qtail) = icell

    do while (qhead <= qtail)
      cur = queue(qhead)
      qhead = qhead + 1
      if (.not. cells(cur)%alive) cycle   ! may already have been processed
      call collapse_cell(cur, queue, qtail, QCAP)
    end do
  end subroutine do_T2

  !------------------------------------------------------------------
  ! Merge ALL of cells(icell)'s vertex columns to a single point (kept
  ! as the first one), globally substitute the others, remove icell,
  ! and enqueue any OTHER cell that this substitution pushed below 3
  ! sides so the caller's cascade loop can clean it up too.
  subroutine collapse_cell(icell, queue, qtail, qcap)
    integer(i4), intent(in)    :: icell
    integer(i4), intent(inout) :: queue(:)
    integer(i4), intent(inout) :: qtail
    integer(i4), intent(in)    :: qcap

    integer(i4) :: n, k, survivor, victim, i
    real(dp) :: sum_api(3), sum_bas(3)

    n = cells(icell)%n
    if (n < 1) then
      cells(icell)%alive = .false.
      cells(icell)%n = 0
      return
    end if

    survivor = cells(icell)%vlist(1)
    sum_api = 0.0_dp; sum_bas = 0.0_dp
    do k = 1, n
      sum_api = sum_api + r_api(:, cells(icell)%vlist(k))
      sum_bas = sum_bas + r_bas(:, cells(icell)%vlist(k))
    end do
    r_api(:, survivor) = sum_api / real(n, dp)
    r_bas(:, survivor) = sum_bas / real(n, dp)

    do k = 2, n
      victim = cells(icell)%vlist(k)
      if (victim == survivor) cycle   ! already-collapsed duplicate within this ring
      do i = 1, n_cell
        if (.not. cells(i)%alive .or. i == icell) cycle
        call ring_replace(i, victim, survivor)
        call collapse_duplicates(i)
        if (cells(i)%n < 3 .and. qtail < qcap) then
          write(*,'(A,I0,A,I0,A)') 'INFO: T2 cascade -- cell ', i, ' pushed to n=', cells(i)%n, ' sides, queued for forced merge.'
          qtail = qtail + 1
          queue(qtail) = i
        end if
      end do
      vert_alive(victim) = .false.
    end do

    cells(icell)%alive = .false.
    cells(icell)%n = 0
  end subroutine collapse_cell

end module mod_T2
