module mod_data
  ! Core data structures for the 3D vertex model.
  !
  ! Geometric model
  ! ---------------
  ! The tissue is a closed shell of polyhedral cells sandwiched between an
  ! outer ("apical") and inner ("basal") sphere-like surface. The two
  ! surfaces share the same combinatorial (topological) network: every
  ! "vertex column" j carries an apical position r_api(:,j) and a basal
  ! position r_bas(:,j) directly "below" it. A cell is a polyhedron whose
  ! apical face is the polygon r_api(:,vlist(1:n)), whose basal face is
  ! the polygon r_bas(:,vlist(1:n)) (same ordering), and whose lateral
  ! faces are the quadrilaterals joining consecutive apical/basal edges.
  ! This is the standard columnar-epithelium simplification of the 3D
  ! vertex model (Honda/Okuda-type) applied to a curved shell.
  !
  ! Every interior vertex column is shared by exactly 3 cells (trivalent),
  ! as guaranteed by the geodesic-icosahedral construction and preserved
  ! by the T1/T2/division rules implemented here.
  use mod_kinds
  implicit none

  integer(i4), parameter :: MAX_SIDES = 16   ! hard cap on polygon valence

  type :: cell_t
    integer(i4) :: n = 0                 ! number of sides (0 = unused slot)
    integer(i4) :: vlist(MAX_SIDES) = 0  ! ordered (CCW, seen from outside) vertex-column ids
    logical     :: alive = .false.
    real(dp)    :: V0 = 0.0_dp           ! target volume
    real(dp)    :: A0 = 0.0_dp           ! target apical area
    real(dp)    :: V_last = 0.0_dp       ! most recent computed volume (set by Force.f90)
    real(dp)    :: A_last = 0.0_dp       ! most recent computed apical area (set by Force.f90)
  end type cell_t

  ! ---- global mutable state ----
  integer(i4) :: n_cell        ! number of cell slots in use (<= size(cells))
  integer(i4) :: n_vert        ! number of vertex-column slots in use (<= NVMAX)
  integer(i4) :: n_cell_cap    ! allocated capacity (cells)
  integer(i4) :: n_vert_cap    ! allocated capacity (vertex columns)

  type(cell_t), allocatable :: cells(:)

  real(dp), allocatable :: r_api(:,:)   ! (3,NVMAX) apical vertex positions
  real(dp), allocatable :: r_bas(:,:)   ! (3,NVMAX) basal vertex positions
  logical,  allocatable :: vert_alive(:)

  real(dp), allocatable :: f_api(:,:)   ! (3,NVMAX) forces on apical vertices
  real(dp), allocatable :: f_bas(:,:)   ! (3,NVMAX) forces on basal vertices

contains

  subroutine allocate_data(ncell_init, nvert_init, growth_factor)
    integer(i4), intent(in) :: ncell_init, nvert_init
    real(dp), intent(in), optional :: growth_factor
    real(dp) :: gf

    gf = 2.0_dp
    if (present(growth_factor)) gf = growth_factor

    n_cell_cap = max(ncell_init, int(ncell_init*gf))
    n_vert_cap = max(nvert_init, int(nvert_init*gf))

    if (allocated(cells))      deallocate(cells)
    if (allocated(r_api))      deallocate(r_api)
    if (allocated(r_bas))      deallocate(r_bas)
    if (allocated(f_api))      deallocate(f_api)
    if (allocated(f_bas))      deallocate(f_bas)
    if (allocated(vert_alive)) deallocate(vert_alive)

    allocate(cells(n_cell_cap))
    allocate(r_api(3, n_vert_cap), r_bas(3, n_vert_cap))
    allocate(f_api(3, n_vert_cap), f_bas(3, n_vert_cap))
    allocate(vert_alive(n_vert_cap))

    vert_alive = .false.
    r_api = 0.0_dp
    r_bas = 0.0_dp
    f_api = 0.0_dp
    f_bas = 0.0_dp

    n_cell = 0
    n_vert = 0
  end subroutine allocate_data

  integer(i4) function count_alive_cells() result(m)
    integer(i4) :: i
    m = 0
    do i = 1, n_cell
      if (cells(i)%alive) m = m + 1
    end do
  end function count_alive_cells

  integer(i4) function count_alive_verts() result(m)
    integer(i4) :: i
    m = 0
    do i = 1, n_vert
      if (vert_alive(i)) m = m + 1
    end do
  end function count_alive_verts

end module mod_data
