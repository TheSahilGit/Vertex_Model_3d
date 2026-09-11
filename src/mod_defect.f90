module mod_defect
  ! Defect-cell protocol (flagged, off by default): seed a random subset
  ! of cells at t=0 with one or more mechanical differences from the
  ! rest of the tissue, then let T4_transition.f90 report separately
  ! whether and in which direction those specific cells get extruded.
  !
  ! Real-world target: the well-documented "apical extrusion of
  ! transformed/abnormal cells" seen in real epithelia and organoids
  ! (epithelial defense against cancer, EDAC) -- a mutant/cancerous cell
  ! surrounded by normal neighbours gets squeezed out of the sheet, and
  ! characteristically toward the APICAL (outer) side.
  !
  ! THREE independent "levers" are provided, each defaulting to a no-op
  ! so any subset can be used alone or combined:
  !
  !   (1) Lambda_line_defect -- this cell's own lateral-tension modulus
  !       (cells(:)%lambda_own, mod_data.f90) instead of Lambda_line.
  !       Tested first, extensively: to actually get CROWDED via this
  !       lever a defect cell needs a HIGHER own lambda than its
  !       neighbours (own lambda_own dominates a cell's own
  !       self-compression -- lower own lambda relaxes to a LARGER
  !       volume, the opposite of "weak and squeezed"), and once
  !       crowded this way its exit direction came out mixed/basal-
  !       leaning across every combination tried (Lambda_line 0.1-0.3,
  !       elevation 0.2-1.2) -- i.e. this lever alone reproduces
  !       preferential removal but NOT the apical bias the real
  !       phenomenon shows.
  !
  !   (2) defect_A0_bas_scale -- scales this cell's OWN basal-area
  !       target (cells(:)%A0_bas) by this extra factor on top of the
  !       usual A0_bas_scale, representing loss of basal identity/
  !       polarity (a well-documented feature of many transformed
  !       epithelial cells, independent of adhesion). Unlike (1), this
  !       is a direct, targeted compressive force on the basal face
  !       specifically (K_A_bas actively pulls Abas toward the smaller
  !       target) while the apical target and Lambda_line are left
  !       exactly as they are for every cell. CONFIRMED (Lambda_line=
  !       Lambda_line_defect=0.3, defect_A0_bas_scale=0.1, threshold
  !       0.975): defect cells are extruded both preferentially (10/16
  !       gone while only 20% of the whole tissue had died) and with a
  !       genuine apical majority (roughly 60% apical across several
  !       thresholds/scales tried) -- the best of the three levers, and
  !       the one that actually reproduces the real phenomenon.
  !
  !   (3) defect_boundary_tension -- an EXTRA lateral tension added only
  !       on edges that connect a defect cell to a non-defect neighbour
  !       (a heterotypic interface), split lambda_line_defect-style as
  !       defect_boundary_tension/2 from each owning cell's own share
  !       (Force.f90). This is the literal "actomyosin purse-string"
  !       mechanism reported for real EDAC. CONFIRMED to behave like
  !       lever (1), as expected from how sparsely defect cells are
  !       normally seeded (most of a defect cell's own edges already
  !       border a normal neighbour, so "boundary-only" tension and
  !       "own-cell" tension act almost identically here): strong
  !       preferential-removal timing (e.g. 9/16 gone while only 5% of
  !       the whole tissue had died, at Lambda_line=Lambda_line_defect=
  !       0.3, defect_boundary_tension=0.5, threshold 0.97), but a
  !       basal-leaning direction (44% apical), not the apical majority
  !       the real phenomenon shows. Included for the comparison, not
  !       as the recommended lever.
  !
  ! Caveat common to all three levers, and to T4 in general at this mesh
  ! size: once T4 activity starts at all, this small (162-cell), tightly
  ! -coupled mesh has consistently cascaded to full tissue loss in every
  ! configuration tried here (a known feedback effect -- extruding a
  ! cell strains its neighbours further -- documented earlier for plain
  ! T4 with no defects at all). The preferential-timing and direction
  ! statistics above are read from the early part of such a run, well
  ! before the bulk tissue starts dying, not from a genuinely steady
  ! "a few cells culled, the rest of the tissue fine forever" state --
  ! no threshold was found that avoids the eventual collapse at this
  ! mesh size. A larger mesh (bigger N_subdivision) or a touch of
  ! thermal noise (kBT>0) are the likeliest ways to get a genuinely
  ! non-collapsing demonstration, if that is wanted.
  !
  ! Mechanism. mesh_init.f90 sets every cell's lambda_own to plain
  ! Lambda_line as it's built, so with defect_enable=.false. this module
  ! is a complete no-op and Force.f90 behaves exactly as it always did.
  ! When enabled, seed_defect_cells -- called once, right after
  ! init_mesh() and right after main.f90 sets the usual t=0 targets,
  ! before the very first real compute_forces of the time loop -- picks
  ! floor(defect_fraction * n_alive) cells uniformly at random (a
  ! partial Fisher-Yates shuffle of the alive-cell index list, using the
  ! same Fortran RNG already seeded by mod_langevin's init_rng), applies
  ! levers (1)-(2) to them directly (lever (3) needs no per-cell state,
  ! it is applied live in Force.f90 keyed off is_defect), and tags them
  ! cells(:)%is_defect = .true. -- a permanent diagnostic marker, never
  ! touched again except by Cell_Division.f90 (both daughters inherit
  ! their parent's lambda_own, A0_bas and is_defect, so a defect cell's
  ! lineage stays defective).
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
      cells(alive_idx(i))%A0_bas     = defect_A0_bas_scale * cells(alive_idx(i))%A0_bas
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
