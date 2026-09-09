function M = read_mesh_meta(fname)
% READ_MESH_META  Read the plain-text data/mesh_meta.txt key/value file
% written once by the Fortran code (mod_io.f90 % write_mesh_meta).
%
%   M = read_mesh_meta()                 % looks for 'data/mesh_meta.txt'
%   M = read_mesh_meta('data/mesh_meta.txt')
%
% Returns a struct with fields MAX_SIDES, R_apical, R_basal, dt,
% it_dumps, n_cell_init, n_vert_init (numeric).

if nargin < 1
    fname = fullfile('data', 'mesh_meta.txt');
end

fid = fopen(fname, 'r');
if fid < 0
    error('read_mesh_meta:open', 'Could not open %s', fname);
end

M = struct();
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    line = strtrim(line);
    if isempty(line) || line(1) == '%'
        continue;
    end
    parts = strsplit(line);
    if numel(parts) < 2
        continue;
    end
    key = parts{1};
    val = str2double(parts{2});
    if isnan(val)
        val = parts{2};
    end
    M.(key) = val;
end
fclose(fid);
end
