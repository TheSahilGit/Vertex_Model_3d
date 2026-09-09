module mod_force
  ! Standard 3D vertex-model energy functional and its analytic gradient
  ! (= minus the mechanical force on every vertex).
  !
  !   E = Sum_cells  (K_V/2)(V_i    - V0_i)^2        volume elasticity
  !     + Sum_cells  (K_A/2)(Aapi_i - A0_i)^2         apical-area elasticity
  !     + Sum_cells  (K_P/2) Papi_i^2                 apical perimeter contractility
  !     + Sum_edges  Lambda_line * S_lateral(edge)    lateral (cell-cell) interfacial tension
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
  ! Each lateral interface is shared by exactly two cells, both of
  ! which loop over "their own" copy of that edge; to avoid double
  ! counting, every cell contributes Lambda_line/2 * S_edge to the
  ! energy (and the corresponding half-weighted gradient) so that the
  ! sum over both owning cells reproduces Lambda_line * S_edge exactly
  ! once.
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  implicit none

contains

  subroutine compute_forces(energy, vol_total, area_total)
    real(dp), intent(out) :: energy
    real(dp), intent(out) :: vol_total, area_total

    integer(i4) :: ic, n, k, kp, j, gvid
    real(dp) :: Aarr(3, MAX_SIDES), Barr(3, MAX_SIDES)
    real(dp) :: ca(3), cb(3)
    real(dp) :: gV_A(3, MAX_SIDES), gV_B(3, MAX_SIDES)
    real(dp) :: gArea_A(3, MAX_SIDES)
    real(dp) :: gPer_A(3, MAX_SIDES)
    real(dp) :: gLat_A(3, MAX_SIDES), gLat_B(3, MAX_SIDES)
    real(dp) :: V, Aapi, Papi, Slat_total
    real(dp) :: vtri, gp(3), gq(3), gr(3)
    real(dp) :: vtri2, gp2(3), gq2(3), gr2(3)
    real(dp) :: atri
    real(dp) :: d(3), edge_len
    real(dp) :: dV, dA

    f_api = 0.0_dp
    f_bas = 0.0_dp
    energy = 0.0_dp
    vol_total = 0.0_dp
    area_total = 0.0_dp

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
      gArea_A(:, 1:n) = 0.0_dp
      gPer_A(:, 1:n) = 0.0_dp
      gLat_A(:, 1:n) = 0.0_dp; gLat_B(:, 1:n) = 0.0_dp
      V = 0.0_dp; Aapi = 0.0_dp; Papi = 0.0_dp; Slat_total = 0.0_dp

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
        Slat_total = Slat_total + vtri + vtri2
        gLat_A(:, k)  = gLat_A(:, k)  + gp + gp2
        gLat_B(:, k)  = gLat_B(:, k)  + gq
        gLat_B(:, kp) = gLat_B(:, kp) + gr + gq2
        gLat_A(:, kp) = gLat_A(:, kp) + gr2
      end do

      dV = V - cells(ic)%V0
      dA = Aapi - cells(ic)%A0

      energy = energy + 0.5_dp * K_V * dV * dV &
                       + 0.5_dp * K_A * dA * dA &
                       + 0.5_dp * K_P * Papi * Papi &
                       + 0.5_dp * Lambda_line * Slat_total
      vol_total  = vol_total  + V
      area_total = area_total + Aapi
      cells(ic)%V_last = V
      cells(ic)%A_last = Aapi

      do k = 1, n
        gvid = cells(ic)%vlist(k)
        f_api(:, gvid) = f_api(:, gvid) - ( K_V * dV * gV_A(:, k) &
                                          + K_A * dA * gArea_A(:, k) &
                                          + K_P * Papi * gPer_A(:, k) &
                                          + 0.5_dp * Lambda_line * gLat_A(:, k) )
        f_bas(:, gvid) = f_bas(:, gvid) - ( K_V * dV * gV_B(:, k) &
                                          + 0.5_dp * Lambda_line * gLat_B(:, k) )
      end do
    end do
  end subroutine compute_forces

end module mod_force
