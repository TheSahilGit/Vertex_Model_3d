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
%       a larger figure) instead of owning a whole figure. Takes
%       precedence over 'Figure' if both are given.
%
%   plot_tissue_3d(..., 'Title', 'Whole tissue')
%       prepend a label above the usual "t = ..." title, so a saved
%       frame/video says which view it is.
%
%   plot_tissue_3d(..., 'ColorMap', 'turbo')
%       any built-in MATLAB colormap name (default 'parula') -- e.g.
%       'turbo', 'jet', 'hot', 'cool', 'copper', 'bone', 'hsv'.
%
%   plot_tissue_3d(..., 'CVal', cval, 'CBLabel', cblabel)
%       use this already-computed per-cell colour vector (and its
%       label) instead of calling tissue_color_values(S,colorby) again
%       -- COLORBY is then ignored for colouring purposes (still used
%       for nothing else, so may be left as whatever). For 'shapefactor'
%       et al., which need cell_lateral_area.m's geometry computation,
%       this lets a whole video (make_tissue_video.m) compute each
%       frame's values ONCE and reuse them, instead of paying for that
%       computation twice per frame (once for the video's global
%       colour-scale scan, once here).
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2 || isempty(colorby)
    colorby = 'nsides';
end
p = inputParser;
p.addParameter('Figure', []);
p.addParameter('Axes', []);
p.addParameter('Title', '');
p.addParameter('CLim', []);
p.addParameter('ColorMap', 'parula');
p.addParameter('CVal', []);
p.addParameter('CBLabel', '');
p.parse(varargin{:});
fig_in  = p.Results.Figure;
ax_in   = p.Results.Axes;
clim_in = p.Results.CLim;
title_label = p.Results.Title;
cmap_name = p.Results.ColorMap;

FONT_SIZE   = 26;  % ~2.5x the MATLAB default (10), per user request
TITLE_SIZE  = 30;

F = cell_faces_matrix(S);
if isempty(p.Results.CVal)
    [cval, cblabel] = tissue_color_values(S, colorby);
else
    cval = p.Results.CVal;
    cblabel = p.Results.CBLabel;
end

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
colormap(ax, cmap_name)

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
    if isempty(title_label)
        ax.Position = [0.03 0.06 0.72 0.88];
    else
        % A 2-line title (panel label + "t = ...") needs more headroom
        % than the single-line case above leaves -- checked directly,
        % it clips the top line otherwise.
        ax.Position = [0.03 0.06 0.72 0.80];
    end
    cb.Units = 'normalized';
    cb.Position = [0.80 0.12 0.045 0.76];
end

% NOW freeze the view for good, via 'axis vis3d' -- but only AFTER the
% axes has its final Position, not before. Freezing earlier locks in
% the zoom/box shape appropriate for the ORIGINAL, not-yet-shrunk axes
% box, which is one failure mode (the sphere and title come out
% mispositioned relative to the smaller box); leaving CameraViewAngle
% on 'auto' is another (MATLAB silently recomputes the zoom on any
% later redraw, seen as the whole sphere abruptly shrinking partway
% through an otherwise perfectly smooth video, no topology event or
% data discontinuity anywhere near it).
%
% 'axis vis3d' (rather than just freezing CameraViewAngleMode by hand,
% as this line used to) is the one that actually matters for anyone
% using this figure interactively: it ALSO freezes PlotBoxAspectRatio,
% not just CameraViewAngle. Leaving PlotBoxAspectRatioMode on 'auto'
% (as plain CameraViewAngleMode='manual' does) lets MATLAB re-fit the
% plot box's shape to the axes region for whatever the CURRENT view
% direction happens to be -- harmless for a fixed, non-interactive
% video frame, but the moment a user drags to rotate the figure
% interactively (e.g. the figure make_tissue_video.m leaves open on its
% last frame), that continuous re-fit is exactly what shows up as the
% whole scene appearing to zoom in and out while it rotates.
axis(ax, 'vis3d');

if isempty(title_label)
    title(ax, sprintf('t = %.4g', S.time), 'FontSize', TITLE_SIZE);
else
    % A 2-line title at the plain single-line TITLE_SIZE overflows the
    % top of the figure (checked directly) -- a smaller size for the
    % 2-line case keeps both lines actually visible.
    title(ax, {title_label, sprintf('t = %.4g', S.time)}, 'FontSize', round(TITLE_SIZE*0.65));
end
end
