module mod_io
  ! Binary snapshot writer. Uses Fortran STREAM access (access='stream')
  ! rather than default sequential unformatted I/O so that NO Fortran
  ! record-length markers are embedded in the file -- the result is a
  ! flat, predictable byte stream that Matlab's fread() can read
  ! directly, with a fixed, documented layout (see Analysis_Code_Matlab
  ! /read_vertex_snapshot.m for the matching reader). All integers are
  ! default kind (4 bytes on every mainstream platform/compiler,
  ! including gfortran) and all reals are real(8) (8-byte double),
  ! matching Matlab's 'int32' and 'double'.
  !
  ! Layout of data/snap_<it8.8>.dat :
  !   int32  it
  !   real8  time
  !   real8  R_apical, R_basal
  !   int32  MAX_SIDES
  !   int32  n_vert          (highest vertex-column slot index in use)
  !   int32  n_cell          (highest cell slot index in use)
  !   -- repeated n_vert times --
  !     int32  alive_flag (0/1)
  !     real8  r_api(3)
  !     real8  r_bas(3)
  !   -- repeated n_cell times --
  !     int32  alive_flag (0/1)
  !     int32  n_sides
  !     int32  vlist(MAX_SIDES)     (1-based vertex-column ids; padded with 0)
  !     real8  V0, A0, V_last, A_last, A0_bas, A_bas_last
  !
  ! Layout of data/diag_<it8.8>.dat (scalar time-series diagnostics,
  ! one small file per dump, same it8.8 naming/cadence as snap_<it8.8>,
  ! so a given timestep's snap_ and diag_ files always pair up):
  !   int32  it
  !   real8  time
  !   real8  energy            (total system energy)
  !   real8  lumen_volume      (volume enclosed by the basal/inner surface)
  !   real8  outer_area        (total apical/outer surface area)
  !   real8  max_force         (largest |force| over every alive vertex,
  !                             apical and basal columns both)
  !   int32  n_cells           (alive cell count)
  !   int32  cumulative_T1     (running total T1 count since t=0)
  !   int32  cumulative_T2     (running total T2 count since t=0)
  !   int32  cumulative_T4            (running total completed T4 extrusions since t=0)
  !   int32  cumulative_T4_apical     (of which classified apical-ward, per T4_direction_deadband)
  !   int32  cumulative_T4_basal      (of which classified basal-ward)
  !   int32  cumulative_T4_ambiguous  (of which too close to call within the deadband)
  use mod_kinds
  use mod_parameters
  use mod_data
  implicit none

contains

  subroutine ensure_data_dir()
    logical :: dir_exists
    inquire(file='data/.', exist=dir_exists)
    if (.not. dir_exists) call execute_command_line('mkdir -p data')
  end subroutine ensure_data_dir

  subroutine write_mesh_meta()
    integer :: iun
    open(newunit=iun, file='data/mesh_meta.txt', status='replace', action='write')
    write(iun,'(A)')            '% mesh_meta.txt -- metadata for Analysis_Code_Matlab (key value)'
    write(iun,'(A,1X,I0)')      'MAX_SIDES', MAX_SIDES
    write(iun,'(A,1X,ES16.8)')  'R_apical', R_apical
    write(iun,'(A,1X,ES16.8)')  'R_basal', R_basal
    write(iun,'(A,1X,ES16.8)')  'dt', dt
    write(iun,'(A,1X,I0)')      'it_dumps', it_dumps
    write(iun,'(A,1X,I0)')      'n_cell_init', n_cell
    write(iun,'(A,1X,I0)')      'n_vert_init', n_vert
    close(iun)
  end subroutine write_mesh_meta

  subroutine write_snapshot(it)
    integer(i4), intent(in) :: it
    character(len=256) :: fname
    integer :: iun, j, i
    integer(i4) :: aflag

    write(fname, '(A,I8.8,A)') 'data/snap_', it, '.dat'
    open(newunit=iun, file=trim(fname), status='replace', access='stream', &
         form='unformatted', action='write')

    write(iun) it
    write(iun) real(it, dp) * dt
    write(iun) R_apical, R_basal
    write(iun) MAX_SIDES
    write(iun) n_vert
    write(iun) n_cell

    do j = 1, n_vert
      aflag = 0
      if (vert_alive(j)) aflag = 1
      write(iun) aflag
      write(iun) r_api(:, j)
      write(iun) r_bas(:, j)
    end do

    do i = 1, n_cell
      aflag = 0
      if (cells(i)%alive) aflag = 1
      write(iun) aflag
      write(iun) cells(i)%n
      write(iun) cells(i)%vlist(1:MAX_SIDES)
      write(iun) cells(i)%V0, cells(i)%A0, cells(i)%V_last, cells(i)%A_last, &
                 cells(i)%A0_bas, cells(i)%A_bas_last
    end do

    close(iun)

    call append_manifest(it, fname)
  end subroutine write_snapshot

  subroutine append_manifest(it, fname)
    integer(i4), intent(in) :: it
    character(len=*), intent(in) :: fname
    integer :: iun
    logical :: exist_flag
    inquire(file='data/dump_list.txt', exist=exist_flag)
    if (exist_flag) then
      open(newunit=iun, file='data/dump_list.txt', status='old', position='append', action='write')
    else
      open(newunit=iun, file='data/dump_list.txt', status='replace', action='write')
      write(iun,'(A)') '% it   filename'
    end if
    write(iun,'(I10,2X,A)') it, trim(fname)
    close(iun)
  end subroutine append_manifest

  subroutine write_diagnostics(it, time, energy, lumen_volume, outer_area, max_force, &
                                n_cells, cumulative_T1, cumulative_T2, &
                                cumulative_T4, cumulative_T4_apical, cumulative_T4_basal, &
                                cumulative_T4_ambiguous)
    integer(i4), intent(in) :: it
    real(dp),    intent(in) :: time, energy, lumen_volume, outer_area, max_force
    integer(i4), intent(in) :: n_cells, cumulative_T1, cumulative_T2
    integer(i4), intent(in) :: cumulative_T4, cumulative_T4_apical
    integer(i4), intent(in) :: cumulative_T4_basal, cumulative_T4_ambiguous
    character(len=256) :: fname
    integer :: iun

    write(fname, '(A,I8.8,A)') 'data/diag_', it, '.dat'
    open(newunit=iun, file=trim(fname), status='replace', access='stream', &
         form='unformatted', action='write')

    write(iun) it
    write(iun) time
    write(iun) energy
    write(iun) lumen_volume
    write(iun) outer_area
    write(iun) max_force
    write(iun) n_cells
    write(iun) cumulative_T1
    write(iun) cumulative_T2
    write(iun) cumulative_T4
    write(iun) cumulative_T4_apical
    write(iun) cumulative_T4_basal
    write(iun) cumulative_T4_ambiguous

    close(iun)

    call append_diag_manifest(it, fname)
  end subroutine write_diagnostics

  subroutine append_diag_manifest(it, fname)
    integer(i4), intent(in) :: it
    character(len=*), intent(in) :: fname
    integer :: iun
    logical :: exist_flag
    inquire(file='data/diag_list.txt', exist=exist_flag)
    if (exist_flag) then
      open(newunit=iun, file='data/diag_list.txt', status='old', position='append', action='write')
    else
      open(newunit=iun, file='data/diag_list.txt', status='replace', action='write')
      write(iun,'(A)') '% it   filename'
    end if
    write(iun,'(I10,2X,A)') it, trim(fname)
    close(iun)
  end subroutine append_diag_manifest

end module mod_io
