function fig = plot_tissue_3d(S, colorby, varargin)
% PLOT_TISSUE_3D  Standard whole-tissue plot: every alive cell's apical
% face, rendered as a shaded polygon on the outer shell.
%
%   plot_tissue_3d(S)              % color by number of polygon sides
%   plot_tissue_3d(S, 'volume')    % color by V_last/V0 - 1 (volume strain)
%   plot_tissue_3d(S, 'area')      % color by Aapi_last/A0 - 1
%   fig = plot_tissue_3d(...)      % returns the figure handle
%
%   plot_tissue_3d(..., 'Figure', fig, 'CLim', [cmin cmax])
%       Draw into an existing figure (cleared first) instead of
%       creating a new one, with an explicit, fixed colour range
%       instead of auto-scaling to this one frame's data. Used by
%       make_tissue_video.m so a whole video has a constant window
%       size (required by VideoWriter) and a constant, comparable
%       colour scale across frames.
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2 || isempty(colorby)
    colorby = 'nsides';
end
p = inputParser;
p.addParameter('Figure', []);
p.addParameter('CLim', []);
p.parse(varargin{:});
fig_in  = p.Results.Figure;
clim_in = p.Results.CLim;

FONT_SIZE   = 26;  % ~2.5x the MATLAB default (10), per user request
TITLE_SIZE  = 30;

idx = find(S.cell_alive);
F = cell_faces_matrix(S);
[cval, cblabel] = tissue_color_values(S, colorby);

if isempty(fig_in)
    fig = figure('Color', 'w');
else
    fig = fig_in;
    figure(fig);
    clf(fig);
end

patch('Faces', F, 'Vertices', S.r_api, ...
      'FaceVertexCData', cval, 'FaceColor', 'flat', ...
      'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);

axis equal vis3d off
view(35, 20)
camlight('headlight'); lighting gouraud; material dull
colormap(parula)

if isempty(clim_in)
    if max(cval) > min(cval)
        clim([min(cval), max(cval)]);
    end
else
    clim(clim_in);
end

cb = colorbar;
cb.Label.String = cblabel;
cb.FontSize = FONT_SIZE;

set(gca, 'FontSize', FONT_SIZE);
% title(sprintf('Tissue at step %d  (t = %.4g),  N_{cell} = %d', ...
%       S.it, S.time, numel(idx)), 'FontSize', TITLE_SIZE); %#ok<UNRCH>
end
