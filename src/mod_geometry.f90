module mod_geometry
  ! Elementary geometric primitives (triangle area, tetrahedron signed
  ! volume) together with their ANALYTIC gradients with respect to each
  ! of the three vertex positions. These are the building blocks used
  ! by Force.f90 to differentiate cell volume / apical area / apical
  ! perimeter / lateral-face area with respect to every vertex.
  use mod_kinds
  implicit none

contains

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  pure function norm3(a) result(m)
    real(dp), intent(in) :: a(3)
    real(dp) :: m
    m = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
  end function norm3

  !---------------------------------------------------------------------
  ! Signed volume of the tetrahedron (origin, p, q, r):
  !   vol = (1/6) p . (q x r)
  ! Used as the elementary term in the divergence-theorem decomposition
  ! of a closed polyhedron's volume: summing this quantity over every
  ! (consistently outward-oriented) triangle of a closed surface gives
  ! the enclosed volume, independent of the reference point (here the
  ! coordinate origin == sphere centre). Gradients are exact (linear
  ! scalar-triple-product derivatives).
  !---------------------------------------------------------------------
  pure subroutine tetra_vol_grad(p, q, r, vol, gp, gq, gr)
    real(dp), intent(in)  :: p(3), q(3), r(3)
    real(dp), intent(out) :: vol, gp(3), gq(3), gr(3)
    vol = dot_product(p, cross(q, r)) / 6.0_dp
    gp  = cross(q, r) / 6.0_dp
    gq  = cross(r, p) / 6.0_dp
    gr  = cross(p, q) / 6.0_dp
  end subroutine tetra_vol_grad

  !---------------------------------------------------------------------
  ! Area of triangle (p,q,r) and its gradient wrt each vertex.
  !   n = (q-p) x (r-p),  A = |n|/2,  nhat = n/|n|
  !   grad_p A = (1/2) nhat x (r - q)
  !   grad_q A = (1/2) nhat x (p - r)
  !   grad_r A = (1/2) nhat x (q - p)
  ! (standard triangle-area-gradient identity; each gradient lies in the
  ! triangle's plane, is perpendicular to the opposite edge, and the
  ! three gradients sum exactly to zero as required by translational
  ! invariance of the area -- checked numerically in sanity_checks.f90)
  !---------------------------------------------------------------------
  pure subroutine tri_area_grad(p, q, r, area, gp, gq, gr, nhat_out)
    real(dp), intent(in)  :: p(3), q(3), r(3)
    real(dp), intent(out) :: area, gp(3), gq(3), gr(3)
    real(dp), intent(out), optional :: nhat_out(3)
    real(dp) :: n(3), nrm, nhat(3)

    n   = cross(q - p, r - p)
    nrm = norm3(n)
    if (nrm > 1.0e-300_dp) then
      nhat = n / nrm
    else
      nhat = 0.0_dp
    end if
    area = 0.5_dp * nrm

    gp = 0.5_dp * cross(nhat, r - q)
    gq = 0.5_dp * cross(nhat, p - r)
    gr = 0.5_dp * cross(nhat, q - p)

    if (present(nhat_out)) nhat_out = nhat
  end subroutine tri_area_grad

  pure function centroid_n(pts, n) result(c)
    integer(i4), intent(in) :: n
    real(dp), intent(in) :: pts(3, n)
    real(dp) :: c(3)
    integer :: k
    c = 0.0_dp
    do k = 1, n
      c = c + pts(:, k)
    end do
    c = c / real(n, dp)
  end function centroid_n

end module mod_geometry
