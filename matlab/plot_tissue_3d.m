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
%   plot_tissue_3d(..., 'Axes', ax, 'CLim', [cmin cmax])
%       Draw into an existing AXES (cleared first, e.g. one subplot of
%       a larger figure) instead of owning a whole figure -- used by
%       make_combined_video.m to place this view alongside the
%       cross-section ring views in one figure. Takes precedence over
%       'Figure' if both are given.
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2 || isempty(colorby)
    colorby = 'nsides';
end
p = inputParser;
p.addParameter('Figure', []);
p.addParameter('Axes', []);
p.addParameter('CLim', []);
p.parse(varargin{:});
fig_in  = p.Results.Figure;
ax_in   = p.Results.Axes;
clim_in = p.Results.CLim;

FONT_SIZE   = 26;  % ~2.5x the MATLAB default (10), per user request
TITLE_SIZE  = 30;

F = cell_faces_matrix(S);
[cval, cblabel] = tissue_color_values(S, colorby);

if ~isempty(ax_in)
    ax = ax_in;
    cla(ax);
    fig = ancestor(ax, 'figure');
else
    if isempty(fig_in)
        fig = figure('Color', 'w');
    else
        fig = fig_in;
        figure(fig);
        clf(fig);
    end
    ax = axes(fig);
end

patch(ax, 'Faces', F, 'Vertices', S.r_api, ...
      'FaceVertexCData', cval, 'FaceColor', 'flat', ...
      'EdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);

axis(ax, 'equal'); axis(ax, 'off');

% Pin the axes to a FIXED reference frame (the apical radius at
% t=0, a constant, rather than this frame's own data extent) instead
% of leaving x/y/zlim on MATLAB's default "auto" (which refits to
% whatever the current vertex coordinates span). The tissue's actual
% apical radius genuinely, continuously shrinks over the run (nothing
% constrains it to stay at R_apical -- that's only the t=0 value); on
% auto limits, MATLAB re-zooms to compensate and keep the sphere
% filling the frame, which both hides real physical shrinkage AND was
% observed to occasionally do that re-zoom in a sudden jump rather
% than smoothly (the sphere abruptly shrinking for the rest of a video
% at some arbitrary point, with no topology event or other data
% discontinuity anywhere near it). Fixing the limits to a constant
% shows the true shrinkage (correct scientifically) and, as a welcome
% side effect, removes the auto-rezoom mechanism that could glitch.
R = S.R_apical * 1.05;
xlim(ax, [-R R]); ylim(ax, [-R R]); zlim(ax, [-R R]);

view(ax, 35, 20)
camlight(ax, 'headlight'); lighting(ax, 'gouraud'); material(ax, 'dull')
colormap(ax, parula)

if isempty(clim_in)
    if max(cval) > min(cval)
        clim(ax, [min(cval), max(cval)]);
    end
else
    clim(ax, clim_in);
end

cb = colorbar(ax);
cb.Label.String = cblabel;
cb.FontSize = FONT_SIZE;
ax.FontSize = FONT_SIZE;
ax.Toolbar.Visible = 'on';  % 'off' to avoid it showing up in exported/captured frames

% colorbar() auto-shrinks the axes to make room for itself based on
% the CURRENT tick labels' width (e.g. "6" vs "-0.07" need different
% widths); when this view owns the whole figure (no Axes given), fix
% both explicitly (normalized units) so the rendered sphere sits in
% the exact same box on every frame regardless of that. When drawing
% into a caller-supplied Axes (a subplot of a bigger figure), leave
% the caller's own subplot layout alone instead.
if isempty(ax_in)
    ax.Units = 'normalized';
    ax.Position = [0.03 0.06 0.72 0.88];
    cb.Units = 'normalized';
    cb.Position = [0.80 0.12 0.045 0.76];
end

% NOW freeze the camera zoom (CameraViewAngle), i.e. what 'axis vis3d'
% would do -- but only AFTER the axes has its final Position, not
% before. With CameraViewAngleMode left on 'auto' (the default),
% MATLAB is free to silently recompute the zoom on any later redraw,
% which happened partway through an otherwise perfectly smooth
% simulation (no topology event, no data discontinuity at that time)
% and showed up as the whole sphere abruptly shrinking for the rest of
% a video. Freezing it here, once, after Position is final, pins the
% zoom for good. (Freezing earlier -- e.g. via 'axis vis3d' before
% Position was set -- locks in the zoom appropriate for the ORIGINAL,
% not-yet-shrunk axes box, which is the OTHER failure mode: the sphere
% and title come out mispositioned relative to the smaller box.)
ax.CameraViewAngleMode = 'manual';

title(ax, sprintf('t = %.4g', S.time), 'FontSize', TITLE_SIZE);
end
