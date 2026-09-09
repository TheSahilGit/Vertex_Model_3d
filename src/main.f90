program main
  ! 3D Vertex Model on a hollow spherical shell.
  ! See README.md for the full model description; para.in
  ! for run parameters; and the header comments of Force.f90,
  ! T1_transition.f90, T2_transition.f90 and Cell_Division.f90 for the
  ! derivation of the energy terms and the topological-rearrangement
  ! rules implemented here.
  use mod_kinds
  use mod_parameters
  use mod_data
  use mod_mesh_init
  use mod_force
  use mod_langevin
  use mod_T1
  use mod_T2
  use mod_division
  use mod_io
  use mod_sanity
  implicit none

  integer(i4) :: it, i, n_t1, n_t2, n_div
  real(dp) :: energy, vol_total, area_total
  logical :: ok
  character(len=64) :: tag

  write(*,'(A)') '===================================================='
  write(*,'(A)') ' 3D Vertex Model on a hollow spherical shell'
  write(*,'(A)') '===================================================='

  call read_parameters('para.in')

  ! ---------------------------------------------------------------
  ! Self-tests: run BEFORE building the real mesh (they use their own
  ! tiny throwaway allocation of the same global arrays; allocate_data
  ! safely re-allocates when init_mesh sets up the real mesh below).
  ! ---------------------------------------------------------------
  call test_geometry_gradients(ok)
  if (.not. ok) then
    write(*,'(A)') 'ABORT: geometry gradient self-test failed -- fix mod_geometry before trusting any results.'
    stop 1
  end if
  call test_cube_volume(ok)
  if (.not. ok) then
    write(*,'(A)') 'ABORT: unit-cube volume/area self-test failed -- fix Force.f90 winding convention.'
    stop 1
  end if

  call init_rng(random_seed)
  call init_mesh()
  call topology_check('after mesh_init')
  call ring_integrity_check('after mesh_init', ok)
  if (.not. ok) stop 1

  ! ---------------------------------------------------------------
  ! Target shapes: V0 = V0_scale * (initial volume), A0 = A0_scale *
  ! (initial apical area), i.e. by default the freshly built geodesic
  ! mesh IS the unstrained reference configuration (zero elastic
  ! energy at t=0; only the lateral-tension term is initially active).
  ! ---------------------------------------------------------------
  call compute_forces(energy, vol_total, area_total)
  do i = 1, n_cell
    if (.not. cells(i)%alive) cycle
    cells(i)%V0 = V0_scale * cells(i)%V_last
    cells(i)%A0 = A0_scale * cells(i)%A_last
  end do
  call global_shell_checks(vol_total, area_total)

  call ensure_data_dir()
  call write_mesh_meta()
  call write_snapshot(0_i4)
  write(*,'(A,I0,A,I0,A)') 'mesh_init: initial dump written (', count_alive_cells(), &
       ' cells, ', count_alive_verts(), ' vertex columns).'

  write(*,'(A)') 'Starting time evolution...'
  do it = 1, n_steps

    call compute_forces(energy, vol_total, area_total)
    call langevin_step()

    if (mod(it, it_topology_check) == 0) then
      call attempt_T1_transitions(n_t1)
      call attempt_T2_transitions(n_t2)
      if (n_t1 > 0 .or. n_t2 > 0) then
        write(tag, '(A,I0)') 'step ', it
        call topology_check(trim(tag))
        call ring_integrity_check(trim(tag), ok)
        if (.not. ok) stop 2
        write(*,'(A,I0,A,I0,A,I0)') '  -> T1 events: ', n_t1, '   T2 events: ', n_t2, '   n_cell(alive)=', count_alive_cells()
      end if
    end if

    if (it_division_check > 0) then
      if (mod(it, it_division_check) == 0) then
        call attempt_cell_divisions(n_div)
        if (n_div > 0) then
          write(tag, '(A,I0)') 'step ', it
          call topology_check(trim(tag))
          call ring_integrity_check(trim(tag), ok)
          if (.not. ok) stop 3
          write(*,'(A,I0,A,I0)') '  -> division events: ', n_div, '   n_cell(alive)=', count_alive_cells()
        end if
      end if
    end if

    if (mod(it, it_dumps) == 0) then
      call write_snapshot(it)
      write(*,'(A,I8,A,ES14.6,A,ES14.6,A,ES14.6,A,I0)') &
           'step ', it, '   E=', energy, '   V_tot=', vol_total, &
           '   A_tot=', area_total, '   n_cell=', count_alive_cells()
    end if

  end do

  call topology_check('final')
  call ring_integrity_check('final', ok)
  write(*,'(A)') 'Done.'

end program main
