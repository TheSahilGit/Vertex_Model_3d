function make_cross_section_video(files, cutaxis, cutvalue, outfile, fps, colorby)
% MAKE_CROSS_SECTION_VIDEO  Step through every snapshot in FILES (a
% cell array of filenames, e.g. from list_snapshots.m/
% select_snapshots.m), rendering the cross-section (ring) view
% (plot_cross_section.m) coloured by COLORBY ('nsides'|'volume'|
% 'area', default 'volume'), and save the result as a video at
% OUTFILE. One figure is used for the whole run (reused every frame,
% left open on the last frame when done so you can look at it -- close
% it yourself when you're done, or before calling this again).
%
% CUTAXIS/CUTVALUE are fixed once, here, for the whole video -- not
% re-chosen per frame -- so every frame is a cut through the same
% plane as the tissue evolves.
%
% FILES may be a single filename (in a 1x1 cell array): the same code
% path runs, writing a 1-frame video -- which is, in effect, just that
% one plot, saved to OUTFILE like any other.
%
%   make_cross_section_video(files, 'y', 0.0, 'videos/cross_y.avi', 8, 'volume')
%
% Uses the 'Motion JPEG AVI' VideoWriter profile (available on every
% platform, unlike 'MPEG-4' which Linux MATLAB does not support).
%
% Frames are captured with print(fig,'-RGBImage'), not getframe -- see
% make_tissue_video.m for why (getframe can win a race against
% MATLAB's own renderer and grab a frame before a layout change has
% actually been composited).

if nargin < 5 || isempty(fps)
    fps = 8;
end
if nargin < 6 || isempty(colorby)
    colorby = 'volume';
end

FRAME_W = 1000;
FRAME_H = 900;

fprintf('make_cross_section_video [%s=%.3g, %s]: reading %d snapshots...\n', ...
        cutaxis, cutvalue, colorby, numel(files));
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

% No explicit 'Visible' here: inherits the caller's current default
% (normally 'on', so one figure window updates frame-by-frame).
fig = figure('Color', 'w', 'Resize', 'off', ...
             'Units', 'pixels', 'Position', [100 100 FRAME_W FRAME_H]);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(Sall)
    plot_cross_section(Sall{k}, cutaxis, cutvalue, 'Figure', fig, ...
                        'ColorBy', colorby, 'CLim', [cmin cmax]);
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
fprintf('make_cross_section_video [%s=%.3g, %s]: wrote %s\n', cutaxis, cutvalue, colorby, outfile);
end
