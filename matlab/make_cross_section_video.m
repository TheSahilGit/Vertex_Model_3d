function make_cross_section_video(files, plane, outfile, fps, colorby, titleLabel, viewAngle, colormapName)
% MAKE_CROSS_SECTION_VIDEO  Step through every snapshot in FILES (a
% cell array of filenames, e.g. from list_snapshots.m/
% select_snapshots.m), rendering the cross-section (ring) view
% (plot_cross_section.m) coloured by COLORBY (see tissue_color_values.m
% for the full list: 'nsides', 'volume'/'volume_abs', 'area'/
% 'apical_area', 'basal_area'/'basal_area_abs', 'lateral_area',
% 'total_area', 'shapefactor'; default 'volume'), and save the result
% as a video at OUTFILE. One figure is used for the whole run (reused
% every frame, left open on the last frame when done so you can look
% at it -- close it yourself when you're done, or before calling this
% again).
%
% PLANE identifies the cutting plane -- see plot_cross_section.m:
%   'x' | 'y' | 'z'          cut through the origin along that axis
%   a {normal, point} struct  an arbitrary plane (e.g. from
%                             ring_planes_for_location.m)
% Fixed once here for the whole video -- not re-chosen per frame -- so
% every frame is a cut through the same plane as the tissue evolves.
%
% FILES may be a single filename (in a 1x1 cell array): the same code
% path runs, writing a 1-frame video -- which is, in effect, just that
% one plot, saved to OUTFILE like any other.
%
%   make_cross_section_video(files, 'y', 'videos/cross_y.avi', 8, 'volume')
%   make_cross_section_video(files, latPlane, 'videos/lat.avi', 8, 'shapefactor', ...
%                             'Latitude ring', [], 'turbo')
%
% TITLELABEL (optional), if given, is prepended above the usual
% "t = ..." title on every frame, so a saved video says which ring it
% is (e.g. 'Latitude ring' / 'Longitude ring').
%
% VIEWANGLE (optional) is an explicit [az el] passed straight through
% to plot_cross_section.m's 'View' option -- e.g. one found by hand
% (rotate the figure, then read it back with [az,el]=view(gca)) --
% used for every frame instead of the default automatic tilt.
%
% COLORMAPNAME (optional, default 'parula') is any built-in MATLAB
% colormap name, e.g. 'turbo', 'jet', 'hot', 'cool', 'copper', 'bone'.
%
% Uses the 'Motion JPEG AVI' VideoWriter profile (available on every
% platform, unlike 'MPEG-4' which Linux MATLAB does not support).
%
% Each frame's own colour VALUES (tissue_color_values.m) are computed
% once, in the same up-front read pass used to find the global colour
% range, and cached (Cvals below) rather than recomputed during the
% render loop -- see make_tissue_video.m's own version of this for why
% (an expensive mode like 'shapefactor', which triangulates every
% cell's lateral walls from raw geometry, would otherwise pay for that
% real computation twice per frame).
%
% Frames are captured with print(fig,'-RGBImage'), not getframe -- see
% make_tissue_video.m for why (getframe can win a race against
% MATLAB's own renderer and grab a frame before a layout change has
% actually been composited).

if nargin < 4 || isempty(fps)
    fps = 8;
end
if nargin < 5 || isempty(colorby)
    colorby = 'volume';
end
if nargin < 6
    titleLabel = '';
end
if nargin < 7
    viewAngle = [];
end
if nargin < 8 || isempty(colormapName)
    colormapName = 'parula';
end

FRAME_W = 1000;
FRAME_H = 900;

fprintf('make_cross_section_video [%s]: reading %d snapshots...\n', colorby, numel(files));
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

% No explicit 'Visible' here: inherits the caller's current default
% (normally 'on', so one figure window updates frame-by-frame).
fig = figure('Color', 'w', 'Resize', 'off', ...
             'Units', 'pixels', 'Position', [100 100 FRAME_W FRAME_H]);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(Sall)
    plot_cross_section(Sall{k}, plane, 'Figure', fig, ...
                        'ColorBy', colorby, 'CLim', [cmin cmax], 'Title', titleLabel, ...
                        'View', viewAngle, 'ColorMap', colormapName, ...
                        'CVal', Cvals{k}, 'CBLabel', cblabel);
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
fprintf('make_cross_section_video [%s]: wrote %s\n', colorby, outfile);
end
