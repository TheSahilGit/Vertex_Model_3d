module mod_T4
  ! T4: crowding-induced live-cell extrusion.
  !
  ! Biological motivation: real epithelia remove cells two distinct ways.
  ! T2 (T2_transition.f90) is the apoptosis-style route -- a cell that has
  ! already been worn down to a triangle by chance/shape collapses to a
  ! single point, and all of its neighbours suddenly meet there. Live-cell
  ! extrusion is different: a cell that is mechanically SQUEEZED by its
  ! neighbours (overcrowding) is pushed out of the sheet while it is still
  ! healthy, and its neighbours close the gap by becoming directly adjacent
  ! to one another -- not by all meeting at one point. This module
  ! implements that second route.
  !
  ! Mechanism (deliberately built from existing, already-validated pieces
  ! -- no new ring-surgery code anywhere in this file):
  !   A cell is "crowded" when its volume has been squeezed well below its
  !   own target, cells(ic)%V_last < V_extrusion_threshold * cells(ic)%V0
  !   -- the exact volume-based mirror of Cell_Division.f90's own growth
  !   criterion (V_last > V_division_threshold * V0), using only fields
  !   already tracked per cell.
  !
  !   A crowded cell with more than 3 sides has ONE of its own edges (its
  !   own currently-shortest one) flipped via T1_transition.f90's
  !   try_T1_edge -- the exact same primitive ordinary T1 neighbour
  !   exchange already uses. A T1 flip's whole effect is that the two
  !   cells NOT currently sharing that edge become neighbours across a new
  !   wall, while the two cells that DID share it each lose a side. Doing
  !   this repeatedly, once per it_topology_check tick, to a crowded
  !   cell's own shrinking boundary is exactly "neighbours merge faces,
  !   one new adjacency at a time" -- the opposite of T2's single
  !   instantaneous point-collapse, and the reason this needs its own
  !   module rather than reusing T2 outright.
  !
  !   Once a crowded cell is down to 3 sides, one more T1 flip is no
  !   longer meaningful (there is no smaller polygon to shrink toward),
  !   so it is handed directly to T2_transition.f90's do_T2 -- which is
  !   already fully generic to any starting side-count (the n==3
  !   restriction lives only in attempt_T2_transitions, never in do_T2
  !   itself) -- to finish the removal the same way an ordinary T2 event
  !   would.
  !
  ! Statelessness: every it_topology_check tick simply re-scans every
  ! alive cell from scratch. A crowded cell gets at most one flip (or one
  ! final do_T2 call) per tick; nothing about "this cell is mid-extrusion"
  ! is remembered anywhere between ticks -- if try_T1_edge declines this
  ! tick (e.g. a momentarily non-trivalent vertex, or a flanking cell
  ! already at MAX_SIDES), the cell is simply reconsidered fresh next
  ! tick, exactly like ordinary T1 already behaves. If a cell's volume
  ! recovers partway through, it just stops being touched -- no cleanup
  ! required. No new per-cell field is added to mod_data.f90 anywhere.
  !
  ! No upper side-count cap is applied to eligibility: since every single
  ! action taken is one already-validated T1 flip (or a final do_T2 only
  ! once n==3), a large polygon can never vanish violently in one step --
  ! it simply costs more ticks, one flip at a time, however large it
  ! starts.
  !
  ! Direction (apical-ward vs basal-ward): classified at the moment of
  ! final removal directly from the ACTUAL mechanical force on the
  ! cell's own two rings -- not from an area-compression proxy. Every
  ! step's compute_forces (Force.f90) already leaves f_api/f_bas holding
  ! the full analytic force (every energy term: volume, apical/basal
  ! area and perimeter, lateral tension) on every vertex column, so no
  ! new force computation is needed here, only a projection and a sum.
  !
  ! Since the shell is centred on the origin, "radial" at a vertex j is
  ! simply r_api(:,j)/|r_api(:,j)| (or r_bas). Summing the radial
  ! component of the real force over a ring's own vertices gives that
  ! ring's net resultant radial force:
  !   Fnet_api = Sum_{j in ring}  f_api(:,j) . rhat(r_api(:,j))
  !   Fnet_bas = Sum_{j in ring}  f_bas(:,j) . rhat(r_bas(:,j))
  ! Whichever ring's net force is more OUTWARD (larger/less negative) is
  ! being pushed out rather than squeezed in -- that is the exit
  ! direction.
  !
  ! This replaces an earlier version of this rule that instead compared
  ! how compressed each face's area is relative to its own target
  ! (ratio_api = A_last/A0 vs ratio_bas = A_bas_last/A0_bas, the larger
  ! ratio -- the less-compressed face -- being taken as the exit side).
  ! That area-ratio rule is physically intuitive and was itself checked
  ! against this code's own analytic forces (tetra_vol_grad/
  ! tri_area_grad, exactly as used in Force.f90) on a synthetic tapered
  ! ("wedge") cell: the lateral-tension term acts like a drawstring
  ! around each ring's own perimeter, and a BIGGER ring has more
  ! perimeter for that drawstring to grip, so it gets cinched in harder
  ! -- like squeezing a party balloon or a toothpaste tube, where
  ! squeezing near the fat end pushes material out the narrow end. The
  ! net-resultant-force rule below was checked against that same
  ! synthetic geometry and agrees with the area-ratio rule in every
  ! case, but it is the more fundamental of the two -- it reads the
  ! mechanical force straight out of the energy functional instead of
  ! inferring direction indirectly from a compressed-area ratio -- so it
  ! is what is actually implemented here.
  !
  ! The deadband compares the two net forces via their FRACTIONAL
  ! difference, (Fnet_api - Fnet_bas) / (|Fnet_api| + |Fnet_bas|), a
  ! dimensionless quantity in [-1,1] -- so T4_direction_deadband keeps
  ! the same meaning and default (0.05) it always had. Near-tie cases
  ! are bucketed as "ambiguous" rather than forcing a noisy call either
  ! way, since f_api/f_bas are instantaneous quantities recomputed every
  ! step.
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  use mod_T1, only: try_T1_edge
  use mod_T2, only: do_T2
  implicit none

  integer(i4), parameter :: T4_DIR_AMBIGUOUS = 0
  integer(i4), parameter :: T4_DIR_APICAL    = 1
  integer(i4), parameter :: T4_DIR_BASAL     = 2

contains

  subroutine attempt_T4_transitions(n_flips, n_extruded, n_extruded_apical, &
                                     n_extruded_basal, n_extruded_ambiguous)
    integer(i4), intent(out) :: n_flips, n_extruded
    integer(i4), intent(out) :: n_extruded_apical, n_extruded_basal, n_extruded_ambiguous
    logical, allocatable :: touched(:)
    integer(i4) :: ic, j, k, dir_code

    n_flips = 0; n_extruded = 0
    n_extruded_apical = 0; n_extruded_basal = 0; n_extruded_ambiguous = 0
    if (.not. T4_enable) return

    allocate(touched(n_vert_cap))
    touched = .false.

    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      if (cells(ic)%V_last >= V_extrusion_threshold * cells(ic)%V0) cycle   ! not crowded

      if (cells(ic)%n == 3) then
        ! Already down to a triangle: one more T1 flip is not meaningful,
        ! finish the removal the same way an ordinary T2 event would.
        dir_code = classify_T4_direction(ic)
        call do_T2(ic)
        n_extruded = n_extruded + 1
        select case (dir_code)
        case (T4_DIR_APICAL)
          n_extruded_apical = n_extruded_apical + 1
        case (T4_DIR_BASAL)
          n_extruded_basal = n_extruded_basal + 1
        case default
          n_extruded_ambiguous = n_extruded_ambiguous + 1
        end select
      else
        ! Still crowded and larger than a triangle: shrink it by one side,
        ! via one ordinary T1 flip of its own shortest edge -- this is the
        ! step where two of its neighbours merge faces and become adjacent.
        call shortest_own_edge(ic, j, k)
        if (j == 0 .or. k == 0) cycle           ! defensive: degenerate ring
        if (touched(j) .or. touched(k)) cycle   ! already mutated earlier this pass
        call try_T1_edge(j, k, touched)
        if (touched(j)) n_flips = n_flips + 1   ! success iff try_T1_edge marked them
      end if
    end do

    deallocate(touched)
  end subroutine attempt_T4_transitions

  !------------------------------------------------------------------
  ! In icell's own ring, find the two (consecutive) vertex-column ids
  ! spanning its currently-shortest edge -- measured the same way
  ! (apical-length) attempt_T1_transitions already measures every edge
  ! in the mesh, so "shortest" means the same thing here as it does
  ! everywhere else in the code.
  subroutine shortest_own_edge(icell, j, k)
    integer(i4), intent(in)  :: icell
    integer(i4), intent(out) :: j, k
    integer(i4) :: n, kk, kp, v1, v2
    real(dp) :: elen, best

    j = 0; k = 0
    n = cells(icell)%n
    if (n < 3) return
    best = huge(1.0_dp)
    do kk = 1, n
      kp = mod(kk, n) + 1
      v1 = cells(icell)%vlist(kk)
      v2 = cells(icell)%vlist(kp)
      elen = norm3(r_api(:, v1) - r_api(:, v2))
      if (elen < best) then
        best = elen
        j = v1; k = v2
      end if
    end do
  end subroutine shortest_own_edge

  !------------------------------------------------------------------
  ! Classify which face a crowded cell is being squeezed out toward,
  ! from the actual net resultant radial force on each of its own two
  ! rings. See the module header for the full reasoning: the ring whose
  ! net radial force is more OUTWARD is the exit direction.
  integer(i4) function classify_T4_direction(icell) result(code)
    integer(i4), intent(in) :: icell
    integer(i4) :: n, kk, j
    real(dp) :: Fnet_api, Fnet_bas, rn, denom, frac

    n = cells(icell)%n
    Fnet_api = 0.0_dp; Fnet_bas = 0.0_dp
    do kk = 1, n
      j = cells(icell)%vlist(kk)
      rn = norm3(r_api(:, j))
      if (rn > 1.0e-300_dp) Fnet_api = Fnet_api + dot_product(f_api(:, j), r_api(:, j)) / rn
      rn = norm3(r_bas(:, j))
      if (rn > 1.0e-300_dp) Fnet_bas = Fnet_bas + dot_product(f_bas(:, j), r_bas(:, j)) / rn
    end do

    denom = abs(Fnet_api) + abs(Fnet_bas)
    if (denom < 1.0e-300_dp) then
      code = T4_DIR_AMBIGUOUS   ! both rings force-free -- no confident call
      return
    end if
    frac = (Fnet_api - Fnet_bas) / denom

    if (frac > T4_direction_deadband) then
      code = T4_DIR_APICAL      ! apical ring's net radial force relatively more outward -> exits apically
    else if (frac < -T4_direction_deadband) then
      code = T4_DIR_BASAL       ! basal ring's net radial force relatively more outward -> exits basally
    else
      code = T4_DIR_AMBIGUOUS   ! within the deadband -- no confident call
    end if
  end function classify_T4_direction

end module mod_T4
