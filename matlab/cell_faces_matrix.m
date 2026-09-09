function [F, idx] = cell_faces_matrix(S, only_alive)
% CELL_FACES_MATRIX  Build a patch()-ready NaN-padded Faces matrix from
% a snapshot struct S (see read_vertex_snapshot.m). Row i of F lists
% the (1-based) vertex ids of cell i's polygon, padded with NaN past
% cell_n(i) sides -- this is exactly the format Matlab's patch() wants
% for a set of polygons with varying vertex counts.
%
%   F = cell_faces_matrix(S)              % alive cells only
%   [F, idx] = cell_faces_matrix(S, false) % all cell slots, plus the
%                                          % row indices kept (idx)

if nargin < 2
    only_alive = true;
end

if only_alive
    idx = find(S.cell_alive);
else
    idx = (1:S.n_cell)';
end

F = double(S.cell_vlist(idx, :));
F(F == 0) = NaN;
end
