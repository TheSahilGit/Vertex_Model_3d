function D = read_diagnostics(fname)
% READ_DIAGNOSTICS  Read one binary scalar-diagnostics snapshot written
% by mod_io.f90 (write_diagnostics), in the same stream/unformatted
% layout style as read_vertex_snapshot.m (no Fortran record markers;
% see the header comment of src/io_module.f90 for the authoritative
% layout). All integers are 4-byte, all reals are 8-byte (double),
% little-endian on every mainstream platform.
%
%   D = read_diagnostics('data/diag_00005000.dat')
%
% Returns a struct:
%   D.it              saved iteration number
%   D.time            simulation time (it * dt)
%   D.energy          total system energy
%   D.lumen_volume    volume enclosed by the basal/inner surface
%   D.outer_area      total apical/outer surface area
%   D.max_force       largest |force| over every alive vertex
%   D.n_cells         alive cell count
%   D.cumulative_T1   running total T1 count since t=0
%   D.cumulative_T2   running total T2 count since t=0

fid = fopen(fname, 'r', 'ieee-le');
if fid < 0
    error('read_diagnostics:open', 'Could not open %s', fname);
end
c = onCleanup(@() fclose(fid));

D.it            = fread(fid, 1, 'int32');
D.time          = fread(fid, 1, 'float64');
D.energy        = fread(fid, 1, 'float64');
D.lumen_volume  = fread(fid, 1, 'float64');
D.outer_area    = fread(fid, 1, 'float64');
D.max_force     = fread(fid, 1, 'float64');
D.n_cells       = fread(fid, 1, 'int32');
D.cumulative_T1 = fread(fid, 1, 'int32');
D.cumulative_T2 = fread(fid, 1, 'int32');
end
