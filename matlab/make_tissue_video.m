function make_tissue_video(files, colorby, outfile, fps)
% MAKE_TISSUE_VIDEO  Step through every snapshot in FILES (a cell array
% of filenames, e.g. from list_snapshots.m/select_snapshots.m),
% rendering the whole-tissue view (plot_tissue_3d.m) coloured by
% COLORBY ('nsides'|'volume'|'area'), and save the result as a video
% at OUTFILE.
%
%   make_tissue_video(files, 'nsides', 'videos/tissue_nsides.avi', 8)
%
% Uses the 'Motion JPEG AVI' VideoWriter profile, which (unlike
% 'MPEG-4') is available on every platform, including Linux. The whole
% snapshot set is read once up front so a single, fixed colour range
% (the true min/max over the whole selection) can be used for every
% frame -- otherwise each frame would auto-scale to its own data and
% the colour scale would flicker/rescale from frame to frame.
%
% The output figure/frame size is pinned explicitly (FRAME_W x
% FRAME_H below) and every captured frame is force-resized to exactly
% that if it comes back any other size. writeVideo errors out (crashes
% the whole run) the moment ANY frame's pixel size differs even
% slightly from the first one written -- which does happen in
% practice (a colorbar tick label gaining/losing a digit can shift the
% rendered axes by a pixel) -- so this resize is a cheap, load-bearing
% safety net, not a cosmetic touch.

if nargin < 4 || isempty(fps)
    fps = 8;
end

FRAME_W = 1000;
FRAME_H = 900;

fprintf('make_tissue_video [%s]: reading %d snapshots...\n', colorby, numel(files));
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

fig = figure('Color', 'w', 'Visible', 'off', 'Resize', 'off', ...
             'Units', 'pixels', 'Position', [100 100 FRAME_W FRAME_H]);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(Sall)
    plot_tissue_3d(Sall{k}, colorby, 'Figure', fig, 'CLim', [cmin cmax]);
    fig.Units = 'pixels';
    fig.Position = [100 100 FRAME_W FRAME_H];  % re-pin every frame, defensively
    drawnow;

    img = getframe(fig).cdata;
    if size(img, 1) ~= FRAME_H || size(img, 2) ~= FRAME_W
        img = imresize(img, [FRAME_H, FRAME_W]);
    end
    writeVideo(v, img);

    if mod(k, 10) == 0 || k == numel(Sall)
        fprintf('  frame %d/%d\n', k, numel(Sall));
    end
end

close(v);
close(fig);
fprintf('make_tissue_video [%s]: wrote %s\n', colorby, outfile);
end
