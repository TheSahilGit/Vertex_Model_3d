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

ax = gca;
axis(ax, 'equal'); axis(ax, 'off');
view(ax, 35, 20)
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
ax.FontSize = FONT_SIZE;

% colorbar() auto-shrinks the axes to make room for itself based on
% the CURRENT tick labels' width (e.g. "6" vs "-0.07" need different
% widths); fix both explicitly (normalized units) so the rendered
% sphere sits in the exact same box on every frame regardless of that.
%
% NOTE: this only reliably takes effect with 'vis3d' NOT set on the
% axes (see above). 'axis vis3d' freezes CameraViewAngle for
% interactive rotation, and doing that BEFORE resizing Position here
% left the two fighting each other unpredictably -- observed as the
% rendered sphere randomly larger/smaller and even the title going
% missing, frame to frame, with no other change and no error. We don't
% need interactive rotation for a scripted, fixed-view video/plot, so
% it's simply left off rather than juggling the two.
ax.Units = 'normalized';
ax.Position = [0.03 0.06 0.72 0.88];
cb.Units = 'normalized';
cb.Position = [0.80 0.12 0.045 0.76];

title(ax, sprintf('t = %.4g', S.time), 'FontSize', TITLE_SIZE);
end
