function make_tissue_video(files, colorby, outfile, fps, colormapName)
% MAKE_TISSUE_VIDEO  Step through every snapshot in FILES (a cell array
% of filenames, e.g. from list_snapshots.m/select_snapshots.m),
% rendering the whole-tissue view (plot_tissue_3d.m) coloured by
% COLORBY (see tissue_color_values.m for the full list: 'nsides',
% 'volume'/'volume_abs', 'area'/'apical_area', 'basal_area'/
% 'basal_area_abs', 'lateral_area', 'total_area', 'shapefactor'), and
% save the result as a video at OUTFILE. One figure is used for the
% whole run (reused every frame, left open on the last frame when done
% so you can look at it -- close it yourself when you're done, or
% before calling this again).
%
% COLORMAPNAME (optional, default 'parula') is any built-in MATLAB
% colormap name, e.g. 'turbo', 'jet', 'hot', 'cool', 'copper', 'bone'.
%
% FILES may be a single filename (in a 1x1 cell array, e.g. from
% select_snapshots.m with a range that matches exactly one saved
% snapshot): the same code path runs, writing a 1-frame video -- which
% is, in effect, just that one plot, saved to OUTFILE like any other.
%
%   make_tissue_video(files, 'nsides', 'videos/tissue_nsides.avi', 8)
%   make_tissue_video(files, 'shapefactor', 'videos/tissue_q.avi', 8, 'turbo')
%
% Uses the 'Motion JPEG AVI' VideoWriter profile, which (unlike
% 'MPEG-4') is available on every platform, including Linux. The whole
% snapshot set is read once up front so a single, fixed colour range
% (the true min/max over the whole selection) can be used for every
% frame -- otherwise each frame would auto-scale to its own data and
% the colour scale would flicker/rescale from frame to frame. Each
% frame's own colour VALUES (tissue_color_values.m) are computed in
% that same up-front pass and cached (Cvals below), then handed
% straight to plot_tissue_3d.m via its 'CVal' option instead of being
% recomputed during the render loop -- for a cheap mode (e.g. 'nsides',
% a stored field lookup) that saves nothing worth mentioning, but for
% 'total_area'/'shapefactor' (which triangulate every cell's lateral
% walls from raw geometry, cell_lateral_area.m) it avoids doing that
% real computation twice per frame.
%
% The output figure/frame size is pinned explicitly (FRAME_W x
% FRAME_H below) and every captured frame is force-resized to exactly
% that if it comes back any other size. writeVideo errors out (crashes
% the whole run) the moment ANY frame's pixel size differs even
% slightly from the first one written -- which does happen in
% practice (a colorbar tick label gaining/losing a digit can shift the
% rendered axes by a pixel) -- so this resize is a cheap, load-bearing
% safety net, not a cosmetic touch.
%
% Frames are captured with print(fig,'-RGBImage'), not getframe: the
% latter grabs whatever is in the figure's frame buffer at that instant
% and is well known to sometimes win the race against MATLAB's own
% renderer -- observed here as the title or the fixed axes/colorbar
% layout (both applied inside plot_tissue_3d.m) missing from a
% frame's capture even though drawnow had already been called.
% print(...,'-RGBImage') renders synchronously (the same pipeline used
% for file export) and does not have that race.

if nargin < 4 || isempty(fps)
    fps = 8;
end
if nargin < 5 || isempty(colormapName)
    colormapName = 'parula';
end

FRAME_W = 1000;
FRAME_H = 900;

fprintf('make_tissue_video [%s]: reading %d snapshots...\n', colorby, numel(files));
Sall  = cell(numel(files), 1);
Cvals = cell(numel(files), 1);
cblabel = '';
cmin = inf; cmax = -inf;
for k = 1:numel(files)
    Sall{k} = read_vertex_snapshot(files{k});
    [Cvals{k}, cblabel] = tissue_color_values(Sall{k}, colorby);
    cmin = min(cmin, min(Cvals{k}));
    cmax = max(cmax, max(Cvals{k}));
end
if cmax <= cmin
    cmax = cmin + 1;
end

% No explicit 'Visible' here: this inherits whatever the caller's
% current default is (normally 'on', so you see one figure window
% update frame-by-frame as the video is built; scripts that want this
% headless can set groot's DefaultFigureVisible to 'off' beforehand).
fig = figure('Color', 'w', 'Resize', 'off', ...
             'Units', 'pixels', 'Position', [100 100 FRAME_W FRAME_H]);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(Sall)
    plot_tissue_3d(Sall{k}, colorby, 'Figure', fig, 'CLim', [cmin cmax], ...
                   'Title', 'Whole tissue', 'ColorMap', colormapName, ...
                   'CVal', Cvals{k}, 'CBLabel', cblabel);
    % Re-pin only the SIZE every frame (defensively) -- never the
    % on-screen location, which is left free for you to drag the
    % window around (e.g. to another monitor) while it renders.
    fig.Units = 'pixels';
    fig.Position(3:4) = [FRAME_W FRAME_H];
    drawnow;

    img = print(fig, '-RGBImage', '-r0');  % '-r0' = match the figure's own
                                            % on-screen pixel size, not a
                                            % fixed default DPI (which
                                            % otherwise returns an image at
                                            % the wrong scale entirely)
    if size(img, 1) ~= FRAME_H || size(img, 2) ~= FRAME_W
        img = imresize(img, [FRAME_H, FRAME_W]);
    end
    writeVideo(v, img);

    if mod(k, 10) == 0 || k == numel(Sall)
        fprintf('  frame %d/%d\n', k, numel(Sall));
    end
end

close(v);
fprintf('make_tissue_video [%s]: wrote %s\n', colorby, outfile);
end
