function make_cross_section_video(files, cutaxis, cutvalue, outfile, fps)
% MAKE_CROSS_SECTION_VIDEO  Step through every snapshot in FILES (a
% cell array of filenames, e.g. from list_snapshots.m/
% select_snapshots.m), rendering the cutaway cross-section view
% (plot_cross_section.m), and save the result as a video at OUTFILE.
% One figure is used for the whole run (reused every frame, so you can
% watch it update as the video is built) and closed again before
% returning -- see make_tissue_video.m for why.
%
% FILES may be a single filename (in a 1x1 cell array): the same code
% path runs, writing a 1-frame video -- which is, in effect, just that
% one plot, saved to OUTFILE like any other.
%
%   make_cross_section_video(files, 'y', 0.0, 'videos/cross_y.avi', 8)
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

FRAME_W = 1000;
FRAME_H = 900;

fprintf('make_cross_section_video [%s=%.3g]: %d frame(s)...\n', cutaxis, cutvalue, numel(files));

% No explicit 'Visible' here: inherits the caller's current default
% (normally 'on', so one figure window updates frame-by-frame).
fig = figure('Color', 'w', 'Resize', 'off', ...
             'Units', 'pixels', 'Position', [100 100 FRAME_W FRAME_H]);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(files)
    S = read_vertex_snapshot(files{k});
    plot_cross_section(S, cutaxis, cutvalue, 'Figure', fig);
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

    if mod(k, 10) == 0 || k == numel(files)
        fprintf('  frame %d/%d\n', k, numel(files));
    end
end

close(v);
close(fig);
fprintf('make_cross_section_video [%s=%.3g]: wrote %s\n', cutaxis, cutvalue, outfile);
end
