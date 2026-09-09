function make_cross_section_video(files, cutaxis, cutvalue, outfile, fps)
% MAKE_CROSS_SECTION_VIDEO  Step through every snapshot in FILES (a
% cell array of filenames, e.g. from list_snapshots.m), rendering the
% cutaway cross-section view (plot_cross_section.m), and save the
% result as a video at OUTFILE.
%
%   make_cross_section_video(files, 'y', 0.0, 'videos/cross_y.avi', 8)
%
% Uses the 'Motion JPEG AVI' VideoWriter profile (available on every
% platform, unlike 'MPEG-4' which Linux MATLAB does not support).

fprintf('make_cross_section_video [%s=%.3g]: %d frames...\n', cutaxis, cutvalue, numel(files));

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 900]);

v = VideoWriter(outfile, 'Motion JPEG AVI');
v.FrameRate = fps;
v.Quality = 90;
open(v);

for k = 1:numel(files)
    S = read_vertex_snapshot(files{k});
    plot_cross_section(S, cutaxis, cutvalue, 'Figure', fig);
    drawnow;
    writeVideo(v, getframe(fig));
    if mod(k, 10) == 0 || k == numel(files)
        fprintf('  frame %d/%d\n', k, numel(files));
    end
end

close(v);
close(fig);
fprintf('make_cross_section_video [%s=%.3g]: wrote %s\n', cutaxis, cutvalue, outfile);
end
