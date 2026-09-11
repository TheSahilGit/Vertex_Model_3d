module mod_defect
  ! Defect-cell protocol (flagged, off by default): seed a random subset
  ! of cells at t=0 with LOWER lateral-tension/adhesion
  ! (Lambda_line_defect) than the rest of the tissue (plain
  ! Lambda_line), then let T4_transition.f90 report separately whether
  ! and in which direction those specific cells get extruded.
  !
  ! Motivation. The T4 direction rule (T4_transition.f90,
  ! classify_T4_direction) was shown to rest on two mechanisms pulling
  ! in OPPOSITE directions: the volume/area/perimeter terms behave like
  ! a "keystone" -- squeezing a wedge-shaped cell pushes it toward its
  ! WIDER face, apical here since R_apical>R_basal -- while the lateral
  ! tension term Lambda acts like a "drawstring" that cinches the WIDER
  ! ring in harder, favouring basal instead. Which one wins depends on
  ! Lambda relative to K_V/K_A: numerically, a single crowded cell in
  ! this shell (K_V=K_A=1, K_P=0.02, K_A_bas=1, K_P_bas=0.1, R_apical=10,
  ! R_basal=8) flips from apical-favouring to basal-favouring somewhere
  ! around Lambda_line~0.6, and every cell tested was apical-favouring,
  ! by a comfortable margin, at Lambda_line=0.3.
  !
  ! A cell with LOWER-than-normal adhesion should sit closer to the
  ! keystone-dominated (apical-favouring) regime than an ordinary
  ! neighbour at the tissue's normal Lambda_line -- this protocol tests
  ! that directly, in situ, inside an otherwise ordinary running
  ! simulation rather than a synthetic single-cell test: seed a few
  ! individually low-adhesion cells inside an otherwise normal tissue,
  ! let the usual dynamics run, and see (a) whether they are actually
  ! the ones that end up crowded and extruded (rather than an ordinary
  ! cell), and (b) which direction they favour, via T4_transition.f90's
  ! new n_extruded_defect* outputs and main.f90's matching
  ! cumulative_T4_defect* counters (data/diag_<it>.dat).
  !
  ! Mechanism. cells(:)%lambda_own (mod_data.f90) already lets every
  ! cell carry its own lateral-tension modulus: mesh_init.f90 sets it to
  ! plain Lambda_line for every cell as its cells are built, so with
  ! defect_enable=.false. this module is a no-op and Force.f90 behaves
  ! exactly as it always did (every cell's own value is identical, so
  ! every shared edge sees plain Lambda_line, as before this feature
  ! existed). When enabled, seed_defect_cells -- called once, right
  ! after init_mesh() and right after main.f90 sets the usual t=0
  ! targets, before the very first real compute_forces of the time
  ! loop -- picks floor(defect_fraction * n_alive) cells uniformly at
  ! random (a partial Fisher-Yates shuffle of the alive-cell index
  ! list, using the same Fortran RNG already seeded by
  ! mod_langevin's init_rng), sets their lambda_own to
  ! Lambda_line_defect, and tags them cells(:)%is_defect = .true. -- a
  ! permanent diagnostic marker, never touched again except by
  ! Cell_Division.f90 (both daughters inherit their parent's lambda_own
  ! and is_defect, so a defect cell's lineage stays low-adhesion).
  use mod_kinds
  use mod_parameters
  use mod_data
  implicit none

contains

  subroutine seed_defect_cells(n_seeded)
    integer(i4), intent(out) :: n_seeded
    integer(i4), allocatable :: alive_idx(:)
    integer(i4) :: n_alive, n_defect, i, j, tmp
    real(dp) :: u

    n_seeded = 0
    if (.not. defect_enable) return

    n_alive = count_alive_cells()
    if (n_alive == 0) return
    allocate(alive_idx(n_alive))
    j = 0
    do i = 1, n_cell
      if (.not. cells(i)%alive) cycle
      j = j + 1
      alive_idx(j) = i
    end do

    n_defect = max(1, nint(defect_fraction * real(n_alive, dp)))
    n_defect = min(n_defect, n_alive)

    ! Partial Fisher-Yates: for i=1..n_defect, swap slot i with a
    ! uniformly random slot in [i, n_alive] -- after n_defect passes,
    ! alive_idx(1:n_defect) is a uniform-random subset of the alive
    ! cells, without replacement.
    do i = 1, n_defect
      call random_number(u)
      j = i + min(int(u * real(n_alive - i + 1, dp)), n_alive - i)
      tmp = alive_idx(i); alive_idx(i) = alive_idx(j); alive_idx(j) = tmp
      cells(alive_idx(i))%lambda_own = Lambda_line_defect
      cells(alive_idx(i))%is_defect  = .true.
    end do

    n_seeded = n_defect
    deallocate(alive_idx)
  end subroutine seed_defect_cells

  !------------------------------------------------------------------
  ! How many seeded defect cells are still alive right now -- an
  ! instantaneous count (not cumulative), written to every
  ! data/diag_<it>.dat dump so a run's defect population over time can
  ! be read straight off the diagnostics without re-deriving it from
  ! the seed list.
  integer(i4) function count_alive_defect_cells() result(m)
    integer(i4) :: i
    m = 0
    do i = 1, n_cell
      if (cells(i)%alive .and. cells(i)%is_defect) m = m + 1
    end do
  end function count_alive_defect_cells

end module mod_defect
