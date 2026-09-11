module mod_force
  ! Standard 3D vertex-model energy functional and its analytic gradient
  ! (= minus the mechanical force on every vertex).
  !
  !   E = Sum_cells  (K_V/2)(V_i        - V0_i)^2         volume elasticity
  !     + Sum_cells  (K_A/2)(Aapi_i     - A0_i)^2         apical-area elasticity
  !     + Sum_cells  (K_P/2) Papi_i^2                     apical perimeter contractility
  !     + Sum_cells  (K_A_bas/2)(Abas_i - A0_bas_i)^2     basal-area elasticity
  !     + Sum_cells  (K_P_bas/2) Pbas_i^2                 basal perimeter contractility
  !     + Sum_edges  Lambda_line * S_lateral(edge)        lateral (cell-cell) interfacial tension
  !       (or, with mod_defect.f90's optional per-cell override active, each
  !       cell uses its OWN cells(ic)%lambda_own for its half-share of every
  !       lateral face it owns -- see that module and the note below)
  !
  ! Each cell is the polyhedron with apical polygon r_api(vlist), basal
  ! polygon r_bas(vlist) and lateral quadrilateral faces joining them.
  ! Volume and areas are computed by triangulating every face (apical
  ! and basal faces are fan-triangulated from their own centroid;
  ! lateral faces are split into 2 triangles) and using the exact
  ! divergence-theorem identity  V = (1/6) Sum_faces (a . (b x c))
  ! over the whole closed, consistently outward-oriented surface
  ! (valid for ANY reference point, here the sphere centre/origin).
  ! Every geometric quantity is differentiated analytically via
  ! mod_geometry (chain rule through the centroid where relevant); see
  ! sanity_checks.f90 for automated checks of both the raw gradients
  ! and the sign convention (winding) used here, against a unit cube.
  !
  ! The basal-area/perimeter terms (K_A_bas, K_P_bas) are the exact
  ! mirror of the apical-area/perimeter ones, just built from the
  ! basal fan/edges instead of the apical ones. Unlike the volume fan,
  ! the basal AREA fan does NOT need a reversed vertex order: area
  ! (|n|/2) and its gradient are both invariant under swapping the
  ! last two triangle corners (swapping flips the sign of nhat, but
  ! that sign cancels out of every gradient formula in
  ! tri_area_grad -- only the signed tetrahedron volume actually cares
  ! about winding). So the same (cb, Barr(kp), Barr(k)) triangle order
  ! already used for the reversed-winding volume fan is reused as-is
  ! for the basal area and perimeter below, with no separate pass.
  !
  ! Each lateral interface is shared by exactly two cells, both of
  ! which loop over "their own" copy of that edge; to avoid double
  ! counting, every cell contributes lambda_own/2 * S_edge to the
  ! energy (and the corresponding half-weighted gradient), where
  ! lambda_own is that CELL's own lateral-tension modulus
  ! (cells(ic)%lambda_own, mod_data.f90) -- equal to Lambda_line for
  ! every cell unless mod_defect.f90 has lowered it for a seeded subset.
  ! So a shared edge's real total tension is (lambda_own_i +
  ! lambda_own_j)/2 * S_edge, the average of its two owning cells' own
  ! values; with no defects seeded this is just Lambda_line * S_edge
  ! exactly as before.
  !
  ! defect_boundary_tension (mod_defect.f90, optional, off/0 by default)
  ! adds a further per-EDGE bonus, ONLY on edges connecting a defect
  ! cell to a non-defect neighbour (a heterotypic interface): each of
  ! the two owning cells adds defect_boundary_tension/2 to its own
  ! share of THAT edge specifically (same split convention as
  ! lambda_own), on top of its own lambda_own/2. This needs to know the
  ! neighbour across each edge; edge_cells (mod_topology.f90) would do
  ! that but costs O(n_cell) PER EDGE, which would make compute_forces
  ! (called every step) effectively O(n_cell^2) -- so when this feature
  ! is active, a cheap vertex-incidence cache (mod_topology.f90's
  ! build_vertex_incidence/neighbor_across_edge) is built once per call
  ! instead, keeping the whole routine at its usual O(n_cell*MAX_SIDES).
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  use mod_topology, only: build_vertex_incidence, neighbor_across_edge
  implicit none

contains

  subroutine compute_forces(energy, vol_total, area_total)
    real(dp), intent(out) :: energy
    real(dp), intent(out) :: vol_total, area_total

    integer(i4) :: ic, n, k, kp, j, gvid, jc
    real(dp) :: Aarr(3, MAX_SIDES), Barr(3, MAX_SIDES)
    real(dp) :: ca(3), cb(3)
    real(dp) :: gV_A(3, MAX_SIDES), gV_B(3, MAX_SIDES)
    real(dp) :: gArea_A(3, MAX_SIDES), gArea_B(3, MAX_SIDES)
    real(dp) :: gPer_A(3, MAX_SIDES), gPer_B(3, MAX_SIDES)
    real(dp) :: gLat_A(3, MAX_SIDES), gLat_B(3, MAX_SIDES)
    real(dp) :: V, Aapi, Papi, Abas, Pbas, Slat_energy
    real(dp) :: vtri, gp(3), gq(3), gr(3)
    real(dp) :: vtri2, gp2(3), gq2(3), gr2(3)
    real(dp) :: atri
    real(dp) :: d(3), edge_len
    real(dp) :: dV, dA, dAbas, w_lat
    logical :: use_boundary_tension
    integer(i4), allocatable :: vic(:,:), vic_count(:)

    f_api = 0.0_dp
    f_bas = 0.0_dp
    energy = 0.0_dp
    vol_total = 0.0_dp
    area_total = 0.0_dp

    use_boundary_tension = defect_enable .and. (abs(defect_boundary_tension) > 0.0_dp)
    if (use_boundary_tension) then
      allocate(vic(3, n_vert_cap), vic_count(n_vert_cap))
      call build_vertex_incidence(vic, vic_count)
    end if

    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      n = cells(ic)%n

      do k = 1, n
        Aarr(:, k) = r_api(:, cells(ic)%vlist(k))
        Barr(:, k) = r_bas(:, cells(ic)%vlist(k))
      end do
      ca = centroid_n(Aarr(:, 1:n), n)
      cb = centroid_n(Barr(:, 1:n), n)

      gV_A(:, 1:n) = 0.0_dp; gV_B(:, 1:n) = 0.0_dp
      gArea_A(:, 1:n) = 0.0_dp; gArea_B(:, 1:n) = 0.0_dp
      gPer_A(:, 1:n) = 0.0_dp; gPer_B(:, 1:n) = 0.0_dp
      gLat_A(:, 1:n) = 0.0_dp; gLat_B(:, 1:n) = 0.0_dp
      V = 0.0_dp; Aapi = 0.0_dp; Papi = 0.0_dp; Abas = 0.0_dp; Pbas = 0.0_dp
      Slat_energy = 0.0_dp

      ! ---- apical cap: fan from ca, volume via tetra(origin,ca,Ak,Akp) ----
      do k = 1, n
        kp = mod(k, n) + 1
        call tetra_vol_grad(ca, Aarr(:, k), Aarr(:, kp), vtri, gp, gq, gr)
        V = V + vtri
        gV_A(:, k)  = gV_A(:, k)  + gq
        gV_A(:, kp) = gV_A(:, kp) + gr
        do j = 1, n
          gV_A(:, j) = gV_A(:, j) + gp / real(n, dp)
        end do

        call tri_area_grad(ca, Aarr(:, k), Aarr(:, kp), atri, gp, gq, gr)
        Aapi = Aapi + atri
        gArea_A(:, k)  = gArea_A(:, k)  + gq
        gArea_A(:, kp) = gArea_A(:, kp) + gr
        do j = 1, n
          gArea_A(:, j) = gArea_A(:, j) + gp / real(n, dp)
        end do

        d = Aarr(:, kp) - Aarr(:, k)
        edge_len = norm3(d)
        Papi = Papi + edge_len
        if (edge_len > 1.0e-14_dp) then
          gPer_A(:, k)  = gPer_A(:, k)  - d / edge_len
          gPer_A(:, kp) = gPer_A(:, kp) + d / edge_len
        end if
      end do

      ! ---- basal cap: fan from cb, REVERSED order -> outward = inward-radial ----
      do k = 1, n
        kp = mod(k, n) + 1
        call tetra_vol_grad(cb, Barr(:, kp), Barr(:, k), vtri, gp, gq, gr)
        V = V + vtri
        gV_B(:, kp) = gV_B(:, kp) + gq
        gV_B(:, k)  = gV_B(:, k)  + gr
        do j = 1, n
          gV_B(:, j) = gV_B(:, j) + gp / real(n, dp)
        end do

        ! basal area: same triangle order as above -- area/its gradient
        ! do not depend on winding (see the header-comment note), so no
        ! separate reversed/non-reversed distinction is needed here.
        call tri_area_grad(cb, Barr(:, kp), Barr(:, k), atri, gp, gq, gr)
        Abas = Abas + atri
        gArea_B(:, kp) = gArea_B(:, kp) + gq
        gArea_B(:, k)  = gArea_B(:, k)  + gr
        do j = 1, n
          gArea_B(:, j) = gArea_B(:, j) + gp / real(n, dp)
        end do

        ! basal perimeter: plain edge length, exactly like the apical one
        d = Barr(:, kp) - Barr(:, k)
        edge_len = norm3(d)
        Pbas = Pbas + edge_len
        if (edge_len > 1.0e-14_dp) then
          gPer_B(:, k)  = gPer_B(:, k)  - d / edge_len
          gPer_B(:, kp) = gPer_B(:, kp) + d / edge_len
        end if
      end do

      ! ---- lateral faces: quad (Ak,Bk,Bkp,Akp) split (Ak,Bk,Bkp)+(Ak,Bkp,Akp) ----
      do k = 1, n
        kp = mod(k, n) + 1

        call tetra_vol_grad(Aarr(:, k), Barr(:, k), Barr(:, kp), vtri, gp, gq, gr)
        V = V + vtri
        gV_A(:, k)  = gV_A(:, k)  + gp
        gV_B(:, k)  = gV_B(:, k)  + gq
        gV_B(:, kp) = gV_B(:, kp) + gr

        call tetra_vol_grad(Aarr(:, k), Barr(:, kp), Aarr(:, kp), vtri2, gp2, gq2, gr2)
        V = V + vtri2
        gV_A(:, k)  = gV_A(:, k)  + gp2
        gV_B(:, kp) = gV_B(:, kp) + gq2
        gV_A(:, kp) = gV_A(:, kp) + gr2

        call tri_area_grad(Aarr(:, k), Barr(:, k), Barr(:, kp), vtri, gp, gq, gr)
        call tri_area_grad(Aarr(:, k), Barr(:, kp), Aarr(:, kp), vtri2, gp2, gq2, gr2)

        ! this face's own effective tension: plain lambda_own, plus a
        ! boundary bonus if this specific edge is a defect/non-defect
        ! heterotypic interface (mod_defect.f90's lever (3))
        w_lat = cells(ic)%lambda_own
        if (use_boundary_tension) then
          jc = neighbor_across_edge(ic, cells(ic)%vlist(k), cells(ic)%vlist(kp), vic, vic_count)
          if (jc > 0) then
            if (cells(ic)%is_defect .neqv. cells(jc)%is_defect) w_lat = w_lat + defect_boundary_tension
          end if
        end if

        Slat_energy = Slat_energy + w_lat * (vtri + vtri2)
        gLat_A(:, k)  = gLat_A(:, k)  + w_lat * (gp + gp2)
        gLat_B(:, k)  = gLat_B(:, k)  + w_lat * gq
        gLat_B(:, kp) = gLat_B(:, kp) + w_lat * (gr + gq2)
        gLat_A(:, kp) = gLat_A(:, kp) + w_lat * gr2
      end do

      dV    = V    - cells(ic)%V0
      dA    = Aapi - cells(ic)%A0
      dAbas = Abas - cells(ic)%A0_bas

      energy = energy + 0.5_dp * K_V * dV * dV &
                       + 0.5_dp * K_A * dA * dA &
                       + 0.5_dp * K_P * Papi * Papi &
                       + 0.5_dp * K_A_bas * dAbas * dAbas &
                       + 0.5_dp * K_P_bas * Pbas * Pbas &
                       + 0.5_dp * Slat_energy
      vol_total  = vol_total  + V
      area_total = area_total + Aapi
      cells(ic)%V_last     = V
      cells(ic)%A_last     = Aapi
      cells(ic)%A_bas_last = Abas

      do k = 1, n
        gvid = cells(ic)%vlist(k)
        f_api(:, gvid) = f_api(:, gvid) - ( K_V * dV * gV_A(:, k) &
                                          + K_A * dA * gArea_A(:, k) &
                                          + K_P * Papi * gPer_A(:, k) &
                                          + 0.5_dp * gLat_A(:, k) )
        f_bas(:, gvid) = f_bas(:, gvid) - ( K_V * dV * gV_B(:, k) &
                                          + K_A_bas * dAbas * gArea_B(:, k) &
                                          + K_P_bas * Pbas * gPer_B(:, k) &
                                          + 0.5_dp * gLat_B(:, k) )
      end do
    end do

    if (use_boundary_tension) deallocate(vic, vic_count)
  end subroutine compute_forces

  !---------------------------------------------------------------------
  ! Volume enclosed by the basal (inner) surface alone -- i.e. the
  ! hollow lumen cavity -- for the "lumen volume" diagnostic.
  !
  ! Uses the same divergence-theorem trick as compute_forces (valid for
  ! any reference point, here the origin/sphere centre, given a closed,
  ! consistently outward-oriented surface), but with the basal fan in
  ! NON-reversed order. Inside compute_forces the basal cap is one face
  ! of a CELL's own closed boundary, so it needs an INWARD (toward the
  ! cell's interior, i.e. away from the lumen) normal -- hence the
  ! reversed fan there. Here the same basal points instead form the
  ! boundary of the LUMEN itself, which needs the opposite convention:
  ! an OUTWARD-from-lumen (away from the origin) normal, i.e. the
  ! non-reversed order used below.
  subroutine compute_lumen_volume(lumen_volume)
    real(dp), intent(out) :: lumen_volume
    integer(i4) :: ic, n, k, kp
    real(dp) :: Barr(3, MAX_SIDES), cb(3)
    real(dp) :: vtri, gp(3), gq(3), gr(3)

    lumen_volume = 0.0_dp
    do ic = 1, n_cell
      if (.not. cells(ic)%alive) cycle
      n = cells(ic)%n
      do k = 1, n
        Barr(:, k) = r_bas(:, cells(ic)%vlist(k))
      end do
      cb = centroid_n(Barr(:, 1:n), n)
      do k = 1, n
        kp = mod(k, n) + 1
        call tetra_vol_grad(cb, Barr(:, k), Barr(:, kp), vtri, gp, gq, gr)
        lumen_volume = lumen_volume + vtri
      end do
    end do
  end subroutine compute_lumen_volume

end module mod_force
