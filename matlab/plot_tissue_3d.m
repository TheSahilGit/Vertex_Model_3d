function fig = plot_tissue_3d(S, colorby)
% PLOT_TISSUE_3D  Standard whole-tissue plot: every alive cell's apical
% face, rendered as a shaded polygon on the outer shell.
%
%   plot_tissue_3d(S)              % color by number of polygon sides
%   plot_tissue_3d(S, 'volume')    % color by V_last/V0 - 1 (volume strain)
%   plot_tissue_3d(S, 'area')      % color by Aapi_last/A0 - 1
%   fig = plot_tissue_3d(...)      % returns the figure handle
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2
    colorby = 'nsides';
end

idx = find(S.cell_alive);
F = cell_faces_matrix(S);

switch lower(colorby)
    case 'nsides'
        cval = S.cell_n(idx);
        cblabel = 'number of sides';
    case 'volume'
        cval = S.cell_Vlast(idx) ./ max(S.cell_V0(idx), eps) - 1;
        cblabel = 'volume strain  (V/V_0 - 1)';
    case 'area'
        cval = S.cell_Alast(idx) ./ max(S.cell_A0(idx), eps) - 1;
        cblabel = 'apical area strain  (A/A_0 - 1)';
    otherwise
        error('plot_tissue_3d:colorby', 'unknown colorby option "%s"', colorby);
end

fig = figure('Color', 'w');
patch('Faces', F, 'Vertices', S.r_api, ...
      'FaceVertexCData', cval, 'FaceColor', 'flat', ...
      'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);

axis equal vis3d off
view(35, 20)
camlight('headlight'); lighting gouraud; material dull
colormap(parula); cb = colorbar; cb.Label.String = cblabel;
title(sprintf('Tissue at step %d  (t = %.4g),  N_{cell} = %d', ...
      S.it, S.time, numel(idx)));
end
