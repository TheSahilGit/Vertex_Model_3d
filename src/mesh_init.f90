module mod_mesh_init
  ! Builds the initial tissue mesh:
  !   1. A geodesic icosahedral triangulation of the unit sphere at
  !      subdivision level L (this plays the role of the "Delaunay"
  !      network -- it is a valid triangulation of the sphere with all
  !      vertices of degree 5 or 6, which is exactly what is needed to
  !      build a well-shaped dual tessellation; it is the standard way
  !      to seed a spherical vertex-model tissue with the physiologically
  !      correct 12 five-fold "defects" and otherwise hexagonal packing).
  !   2. The dual tessellation of that triangulation: one cell per
  !      triangulation VERTEX, one vertex-column per triangulation FACE.
  !      A dual vertex (for triangle (v1,v2,v3) on the unit sphere) is
  !      placed along n = normalize((v2-v1)x(v3-v1)) [sign fixed to the
  !      face's own hemisphere], because this is exactly the point that
  !      is equidistant (in chordal distance) from v1,v2,v3: for unit
  !      vectors, |p-vi|^2 = 2 - 2 p.vi, so p.v1=p.v2=p.v3 requires p to
  !      be parallel to (v2-v1)x(v3-v1). The apical/basal shells are
  !      this same direction scaled by R_apical / R_basal.
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_geometry
  implicit none
  private
  public :: init_mesh

contains

  subroutine init_mesh()
    integer(i4) :: L, NV0, NF0
    real(dp), allocatable :: vtx(:,:)
    integer(i4), allocatable :: faces(:,:)
    integer(i4) :: nv_cur, nf_cur
    integer(i4), allocatable :: edge_mid(:,:)
    integer(i4) :: i, f, v1, v2, v3
    integer(i4), allocatable :: inc_faces(:,:), deg(:)
    integer(i4) :: max_deg_seen

    L = N_subdivision
    NV0 = 10 * 4**L + 2
    NF0 = 20 * 4**L

    write(*,'(A,I0,A,I0,A,I0)') 'mesh_init: L=', L, '  NV0(=n_cell)=', NV0, '  NF0(=n_vert)=', NF0

    allocate(vtx(3, NV0))
    allocate(faces(3, NF0))
    allocate(edge_mid(NV0, NV0))
    edge_mid = 0

    call build_icosahedron(vtx, faces, nv_cur, nf_cur)
    if (nv_cur /= 12 .or. nf_cur /= 20) then
      write(*,*) 'ERROR: base icosahedron malformed', nv_cur, nf_cur
      stop 1
    end if

    do i = 1, L
      call subdivide_once(vtx, faces, nv_cur, nf_cur, edge_mid, NV0, NF0)
    end do

    ! ---- sanity check: exact vertex/face counts for geodesic icosphere ----
    if (nv_cur /= NV0 .or. nf_cur /= NF0) then
      write(*,*) 'ERROR: subdivision produced wrong counts:', nv_cur, NV0, nf_cur, NF0
      stop 1
    end if

    call fix_face_orientation(vtx, faces, nf_cur)

    ! ---- allocate the working (dual / tissue) arrays ----
    call allocate_data(NV0, NF0, growth_factor=capacity_growth_factor)
    n_cell = NV0
    n_vert = NF0

    ! ---- dual vertex positions (one per original face) ----
    do f = 1, nf_cur
      v1 = faces(1, f); v2 = faces(2, f); v3 = faces(3, f)
      call dual_vertex_position(vtx(:, v1), vtx(:, v2), vtx(:, v3), r_api(:, f), r_bas(:, f))
      vert_alive(f) = .true.
    end do

    ! ---- incident-face lists per original vertex (max degree 6) ----
    allocate(inc_faces(6, NV0), deg(NV0))
    deg = 0
    max_deg_seen = 0
    do f = 1, nf_cur
      do i = 1, 3
        v1 = faces(i, f)
        deg(v1) = deg(v1) + 1
        if (deg(v1) > 6) then
          write(*,*) 'ERROR: vertex degree exceeds 6 at original vertex', v1
          stop 1
        end if
        inc_faces(deg(v1), v1) = f
        max_deg_seen = max(max_deg_seen, deg(v1))
      end do
    end do

    ! ---- build each cell = angularly-sorted ring of incident faces ----
    do v1 = 1, NV0
      call build_cell_ring(v1, vtx(:, v1), inc_faces(:, v1), deg(v1))
    end do

    call sanity_check_degrees(deg, NV0)

    write(*,'(A)') 'mesh_init: done.'
  end subroutine init_mesh

  !=====================================================================
  subroutine build_icosahedron(vtx, faces, nv, nf)
    real(dp), intent(inout) :: vtx(:,:)
    integer(i4), intent(inout) :: faces(:,:)
    integer(i4), intent(out) :: nv, nf
    real(dp) :: phi
    real(dp) :: base(3,12)
    integer(i4) :: fl(3,20)
    integer(i4) :: i

    phi = (1.0_dp + sqrt(5.0_dp)) / 2.0_dp

    base(:,1)  = [-1.0_dp,  phi, 0.0_dp]
    base(:,2)  = [ 1.0_dp,  phi, 0.0_dp]
    base(:,3)  = [-1.0_dp, -phi, 0.0_dp]
    base(:,4)  = [ 1.0_dp, -phi, 0.0_dp]
    base(:,5)  = [ 0.0_dp, -1.0_dp,  phi]
    base(:,6)  = [ 0.0_dp,  1.0_dp,  phi]
    base(:,7)  = [ 0.0_dp, -1.0_dp, -phi]
    base(:,8)  = [ 0.0_dp,  1.0_dp, -phi]
    base(:,9)  = [ phi, 0.0_dp, -1.0_dp]
    base(:,10) = [ phi, 0.0_dp,  1.0_dp]
    base(:,11) = [-phi, 0.0_dp, -1.0_dp]
    base(:,12) = [-phi, 0.0_dp,  1.0_dp]

    do i = 1, 12
      vtx(:, i) = base(:, i) / norm3(base(:, i))
    end do

    fl(:,1)  = [1,12,6];  fl(:,2)  = [1,6,2];   fl(:,3)  = [1,2,8]
    fl(:,4)  = [1,8,11];  fl(:,5)  = [1,11,12]
    fl(:,6)  = [2,6,10];  fl(:,7)  = [6,12,5];  fl(:,8)  = [12,11,3]
    fl(:,9)  = [11,8,7];  fl(:,10) = [8,2,9]
    fl(:,11) = [4,10,5];  fl(:,12) = [4,5,3];   fl(:,13) = [4,3,7]
    fl(:,14) = [4,7,9];   fl(:,15) = [4,9,10]
    fl(:,16) = [5,10,6];  fl(:,17) = [3,5,12];  fl(:,18) = [7,3,11]
    fl(:,19) = [9,7,8];   fl(:,20) = [10,9,2]

    do i = 1, 20
      faces(:, i) = fl(:, i)
    end do

    nv = 12
    nf = 20
  end subroutine build_icosahedron

  !=====================================================================
  subroutine subdivide_once(vtx, faces, nv_cur, nf_cur, edge_mid, NV0, NF0)
    real(dp), intent(inout)    :: vtx(:,:)
    integer(i4), intent(inout) :: faces(:,:)
    integer(i4), intent(inout) :: nv_cur, nf_cur
    integer(i4), intent(inout) :: edge_mid(:,:)
    integer(i4), intent(in)    :: NV0, NF0

    integer(i4), allocatable :: new_faces(:,:)
    integer(i4) :: nf_new, f, v1, v2, v3, a, b, c

    allocate(new_faces(3, nf_cur * 4))
    nf_new = 0

    do f = 1, nf_cur
      v1 = faces(1, f); v2 = faces(2, f); v3 = faces(3, f)
      call get_or_create_midpoint(vtx, edge_mid, nv_cur, v1, v2, a, NV0)
      call get_or_create_midpoint(vtx, edge_mid, nv_cur, v2, v3, b, NV0)
      call get_or_create_midpoint(vtx, edge_mid, nv_cur, v3, v1, c, NV0)

      nf_new = nf_new + 1; new_faces(:, nf_new) = [v1, a, c]
      nf_new = nf_new + 1; new_faces(:, nf_new) = [v2, b, a]
      nf_new = nf_new + 1; new_faces(:, nf_new) = [v3, c, b]
      nf_new = nf_new + 1; new_faces(:, nf_new) = [a, b, c]
    end do

    if (nf_new > NF0) then
      write(*,*) 'ERROR: face overflow during subdivision', nf_new, NF0
      stop 1
    end if

    faces(:, 1:nf_new) = new_faces(:, 1:nf_new)
    nf_cur = nf_new
    deallocate(new_faces)
  end subroutine subdivide_once

  !=====================================================================
  subroutine get_or_create_midpoint(vtx, edge_mid, nv_cur, va, vb, id, NV0)
    real(dp), intent(inout)    :: vtx(:,:)
    integer(i4), intent(inout) :: edge_mid(:,:)
    integer(i4), intent(inout) :: nv_cur
    integer(i4), intent(in)    :: va, vb, NV0
    integer(i4), intent(out)   :: id
    integer(i4) :: i, j
    real(dp) :: mid(3)

    i = min(va, vb); j = max(va, vb)
    if (edge_mid(i, j) /= 0) then
      id = edge_mid(i, j)
      return
    end if

    mid = 0.5_dp * (vtx(:, va) + vtx(:, vb))
    mid = mid / norm3(mid)

    nv_cur = nv_cur + 1
    if (nv_cur > NV0) then
      write(*,*) 'ERROR: vertex overflow during subdivision', nv_cur, NV0
      stop 1
    end if
    vtx(:, nv_cur) = mid
    edge_mid(i, j) = nv_cur
    id = nv_cur
  end subroutine get_or_create_midpoint

  !=====================================================================
  subroutine fix_face_orientation(vtx, faces, nf)
    real(dp), intent(in)       :: vtx(:,:)
    integer(i4), intent(inout) :: faces(:,:)
    integer(i4), intent(in)    :: nf
    integer(i4) :: f, v1, v2, v3, tmp
    real(dp) :: n(3), c(3)
    integer(i4) :: n_flipped

    n_flipped = 0
    do f = 1, nf
      v1 = faces(1, f); v2 = faces(2, f); v3 = faces(3, f)
      n = cross(vtx(:, v2) - vtx(:, v1), vtx(:, v3) - vtx(:, v1))
      c = (vtx(:, v1) + vtx(:, v2) + vtx(:, v3)) / 3.0_dp
      if (dot_product(n, c) < 0.0_dp) then
        tmp = faces(2, f); faces(2, f) = faces(3, f); faces(3, f) = tmp
        n_flipped = n_flipped + 1
      end if
    end do
    if (n_flipped > 0) write(*,'(A,I0,A)') 'mesh_init: fixed orientation of ', n_flipped, ' faces'
  end subroutine fix_face_orientation

  !=====================================================================
  subroutine dual_vertex_position(v1, v2, v3, p_api, p_bas)
    real(dp), intent(in)  :: v1(3), v2(3), v3(3)
    real(dp), intent(out) :: p_api(3), p_bas(3)
    real(dp) :: n(3), c(3)
    n = cross(v2 - v1, v3 - v1)
    n = n / norm3(n)
    c = v1 + v2 + v3
    if (dot_product(n, c) < 0.0_dp) n = -n
    p_api = n * R_apical
    p_bas = n * R_basal
  end subroutine dual_vertex_position

  !=====================================================================
  ! Build cell(icell)%vlist as the angularly-sorted cyclic list of the
  ! deg incident face-ids (dual vertex ids), then orient it so the fan
  ! normal points outward (matches r_api direction), i.e. CCW as seen
  ! from outside the shell.
  subroutine build_cell_ring(icell, vhat, inc, deg)
    integer(i4), intent(in) :: icell
    real(dp), intent(in)    :: vhat(3)
    integer(i4), intent(in) :: inc(:)
    integer(i4), intent(in) :: deg
    real(dp) :: e1(3), e2(3), ref(3)
    real(dp) :: ang(6), proj(3)
    integer(i4) :: order(6), i, j, tmpi
    real(dp) :: tmpr
    real(dp) :: n(3)

    if (abs(vhat(3)) < 0.9_dp) then
      ref = [0.0_dp, 0.0_dp, 1.0_dp]
    else
      ref = [1.0_dp, 0.0_dp, 0.0_dp]
    end if
    e1 = cross(vhat, ref); e1 = e1 / norm3(e1)
    e2 = cross(vhat, e1)

    do i = 1, deg
      proj = r_api(:, inc(i)) - dot_product(r_api(:, inc(i)), vhat) * vhat
      ang(i) = atan2(dot_product(proj, e2), dot_product(proj, e1))
      order(i) = inc(i)
    end do

    ! simple insertion sort by angle (deg <= 6)
    do i = 2, deg
      tmpr = ang(i); tmpi = order(i)
      j = i - 1
      do while (j >= 1)
        if (ang(j) <= tmpr) exit
        ang(j+1) = ang(j); order(j+1) = order(j)
        j = j - 1
      end do
      ang(j+1) = tmpr; order(j+1) = tmpi
    end do

    cells(icell)%n = deg
    cells(icell)%vlist = 0
    cells(icell)%vlist(1:deg) = order(1:deg)
    cells(icell)%alive = .true.
    cells(icell)%lambda_own = Lambda_line   ! plain default; mod_defect.f90 may lower this
                                             ! for a random subset, after init_mesh() returns
    cells(icell)%is_defect = .false.

    ! orient outward: test fan normal of first 3 ring vertices
    n = cross(r_api(:, cells(icell)%vlist(2)) - r_api(:, cells(icell)%vlist(1)), &
              r_api(:, cells(icell)%vlist(3)) - r_api(:, cells(icell)%vlist(1)))
    if (dot_product(n, vhat) < 0.0_dp) then
      call reverse_ring(icell)
    end if
  end subroutine build_cell_ring

  subroutine reverse_ring(icell)
    integer(i4), intent(in) :: icell
    integer(i4) :: n, tmp(MAX_SIDES), i
    n = cells(icell)%n
    tmp(1:n) = cells(icell)%vlist(1:n)
    do i = 1, n
      cells(icell)%vlist(i) = tmp(n - i + 1)
    end do
  end subroutine reverse_ring

  !=====================================================================
  subroutine sanity_check_degrees(deg, NV0)
    integer(i4), intent(in) :: deg(:), NV0
    integer(i4) :: i, n5, n6, nother
    n5 = 0; n6 = 0; nother = 0
    do i = 1, NV0
      select case (deg(i))
      case (5); n5 = n5 + 1
      case (6); n6 = n6 + 1
      case default; nother = nother + 1
      end select
    end do
    write(*,'(A,I0,A,I0,A,I0)') 'mesh_init: degree-5 cells=', n5, '  degree-6 cells=', n6, '  other=', nother
    if (n5 /= 12 .or. nother /= 0) then
      write(*,*) 'WARNING: expected exactly 12 pentagonal cells and 0 others (geodesic icosphere property).'
    end if
  end subroutine sanity_check_degrees

end module mod_mesh_init
