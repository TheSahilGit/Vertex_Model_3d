function S = read_vertex_snapshot(fname)
% READ_VERTEX_SNAPSHOT  Read one binary tissue snapshot written by
% mod_io.f90 (write_snapshot), in stream/unformatted layout (no Fortran
% record markers -- see the header comment of src/io_module.f90 for the
% authoritative layout description). All integers are 4-byte, all reals
% are 8-byte (double), little-endian on every mainstream platform.
%
%   S = read_vertex_snapshot('data/snap_00000100.dat')
%
% Returns a struct:
%   S.it, S.time, S.R_apical, S.R_basal, S.MAX_SIDES, S.n_vert, S.n_cell
%   S.vert_alive (n_vert x 1 logical)
%   S.r_api, S.r_bas            (n_vert x 3, apical/basal positions)
%   S.cell_alive (n_cell x 1 logical)
%   S.cell_n                    (n_cell x 1, number of sides)
%   S.cell_vlist                (n_cell x MAX_SIDES, 1-based ids into
%                                 r_api/r_bas; unused slots are 0)
%   S.cell_V0, S.cell_A0        (n_cell x 1, target volume/apical area)
%   S.cell_Vlast, S.cell_Alast  (n_cell x 1, last computed volume/area)
%   S.cell_A0bas, S.cell_Abaslast (n_cell x 1, target/last basal area)

fid = fopen(fname, 'r', 'ieee-le');
if fid < 0
    error('read_vertex_snapshot:open', 'Could not open %s', fname);
end
c = onCleanup(@() fclose(fid));

S.it        = fread(fid, 1, 'int32');
S.time      = fread(fid, 1, 'float64');
rr          = fread(fid, 2, 'float64');
S.R_apical  = rr(1);
S.R_basal   = rr(2);
S.MAX_SIDES = fread(fid, 1, 'int32');
S.n_vert    = fread(fid, 1, 'int32');
S.n_cell    = fread(fid, 1, 'int32');

S.vert_alive = false(S.n_vert, 1);
S.r_api = zeros(S.n_vert, 3);
S.r_bas = zeros(S.n_vert, 3);
for j = 1:S.n_vert
    a = fread(fid, 1, 'int32');
    S.vert_alive(j) = (a ~= 0);
    S.r_api(j, :) = fread(fid, 3, 'float64')';
    S.r_bas(j, :) = fread(fid, 3, 'float64')';
end

S.cell_alive   = false(S.n_cell, 1);
S.cell_n       = zeros(S.n_cell, 1);
S.cell_vlist   = zeros(S.n_cell, S.MAX_SIDES);
S.cell_V0      = zeros(S.n_cell, 1);
S.cell_A0      = zeros(S.n_cell, 1);
S.cell_Vlast   = zeros(S.n_cell, 1);
S.cell_Alast   = zeros(S.n_cell, 1);
S.cell_A0bas    = zeros(S.n_cell, 1);
S.cell_Abaslast = zeros(S.n_cell, 1);
for i = 1:S.n_cell
    a = fread(fid, 1, 'int32');
    S.cell_alive(i) = (a ~= 0);
    S.cell_n(i) = fread(fid, 1, 'int32');
    S.cell_vlist(i, :) = fread(fid, S.MAX_SIDES, 'int32')';
    vals = fread(fid, 6, 'float64');
    S.cell_V0(i)    = vals(1);
    S.cell_A0(i)    = vals(2);
    S.cell_Vlast(i) = vals(3);
    S.cell_Alast(i) = vals(4);
    S.cell_A0bas(i)    = vals(5);
    S.cell_Abaslast(i) = vals(6);
end
end
