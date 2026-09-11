function make_combined_video(files, ref, colorby, outfile, fps)
% MAKE_COMBINED_VIDEO  Step through every snapshot in FILES (a cell
% array of filenames, e.g. from list_snapshots.m/select_snapshots.m),
% rendering THREE views of the same frame side by side in one figure --
% the whole tissue (plot_tissue_3d.m), the LATITUDE ring through a
% given cell/location (plot_cross_section.m), and the LONGITUDE ring
% through that same location -- and saving the result as a single
% video at OUTFILE. One figure (with 3 fixed subplot axes, reused every
% frame) is used for the whole run, left open on the last frame when
% done.
%
% REF identifies the reference cell/location the two rings are cut
% through -- see ring_planes_for_location.m:
%   ref = 7           % a cell index (that cell's own apical centroid)
%   ref = [x y z]      % an explicit point in space
% It is resolved into two FIXED cutting planes ONCE, from the FIRST
% frame in FILES, and those same two planes are reused unchanged for
% every later frame (exactly like the old cut_axis/cut_value were
% fixed for a whole video) -- the rings are not re-centred frame to
% frame even if the tissue moves or that cell's own shape drifts.
%
% COLORBY is 'nsides'|'volume'|'area' (default 'volume'), shared by all
% three panels with one common, fixed colour scale (the true min/max
% over the whole selection) so all three are directly comparable and
% the scale does not flicker frame to frame.
%
% FILES may be a single filename (in a 1x1 cell array): the same code
% path runs, writing a 1-frame video -- which is, in effect, just that
% one 3-panel figure, saved to OUTFILE like any other.
%
%   make_combined_video(files, 7, 'volume', 'videos/combined.avi', 8)
%
% Uses the 'Motion JPEG AVI' VideoWriter profile (available on every
% platform, unlike 'MPEG-4' which Linux MATLAB does not support).
%
% Frames are captured with print(fig,'-RGBImage'), not getframe -- see
% make_tissue_video.m for why (getframe can win a race against
% MATLAB's own renderer and grab a frame before a layout change has
% actually been composited).

if nargin < 3 || isempty(colorby)
    colorby = 'volume';
end
if nargin < 5 || isempty(fps)
    fps = 8;
end

FRAME_W = 2700;
FRAME_H = 900;

fprintf('make_combined_video [%s]: reading %d snapshots...\n', colorby, numel(files));
Sall = cell(numel(files), 1);
cmin = inf; cmax = -inf;
for k = 1:numel(files)
    Sall{k} = read_vertex_snapshot(files{k});
    cval = tissue_color_values(Sall{k}, colorby);
    cmin = min(cmin, min(cval));
    cmax = max(cmax, max(cval));
end
if cmax <= cmin
    cmax = cmin + 1;
end
clims = [cmin cmax];

[latPlane, lonPlane, refPoint] = ring_planes_for_location(Sall{1}, ref);
fprintf('  reference point = [%.4g %.4g %.4g]\n', refPoint);

% No explicit 'Visible' here: inherits the caller's current default
% (normally 'on', so one figure window updates frame-by-frame).
fig = figure('Color', 'w', 'Resize', 'off', ...
             'Units', 'pixels', 'Position', [100 100 FRAME_W FRAME_H]);
tl = tiledlayout(fig, 1, 3, 'Padding', 'loose', 'TileSpacing', 'compact');
PANEL_TITLE_SIZE = 20;  % the individual plot functions already set their
                        % own (larger) single-line "t = ..." title; the
                        % 2-line {panel label; t=...} title set below,
                        % every frame, needs a smaller size or it overflows
                        % the tile's own headroom at the top of the figure
ax_tissue = nexttile(tl, 1);
ax_lat    = nexttile(tl, 2);
ax_lon    = nexttile(tl, 3);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(Sall)
    plot_tissue_3d(Sall{k}, colorby, 'Axes', ax_tissue, 'CLim', clims);
    title(ax_tissue, {'Whole tissue', sprintf('t = %.4g', Sall{k}.time)}, 'FontSize', PANEL_TITLE_SIZE);

    plot_cross_section(Sall{k}, latPlane, 'Axes', ax_lat, 'ColorBy', colorby, 'CLim', clims);
    title(ax_lat, {'Latitude ring', sprintf('t = %.4g', Sall{k}.time)}, 'FontSize', PANEL_TITLE_SIZE);

    plot_cross_section(Sall{k}, lonPlane, 'Axes', ax_lon, 'ColorBy', colorby, 'CLim', clims);
    title(ax_lon, {'Longitude ring', sprintf('t = %.4g', Sall{k}.time)}, 'FontSize', PANEL_TITLE_SIZE);

    % Re-pin only the SIZE every frame (defensively) -- never the
    % on-screen location, which is left free for you to drag the
    % window around (e.g. to another monitor) while it renders.
    fig.Units = 'pixels';
    fig.Position(3:4) = [FRAME_W FRAME_H];
    drawnow;

    img = print(fig, '-RGBImage', '-r0');  % '-r0' = match the figure's own
                                            % on-screen pixel size, not a
                                            % fixed default DPI
    if size(img, 1) ~= FRAME_H || size(img, 2) ~= FRAME_W
        img = imresize(img, [FRAME_H, FRAME_W]);
    end
    writeVideo(v, img);

    if mod(k, 10) == 0 || k == numel(Sall)
        fprintf('  frame %d/%d\n', k, numel(Sall));
    end
end

close(v);
fprintf('make_combined_video [%s]: wrote %s\n', colorby, outfile);
end
