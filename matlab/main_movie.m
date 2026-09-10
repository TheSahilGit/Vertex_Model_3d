% main_movie.m
% Driver: reads the tissue snapshots written by the Fortran code into
% data/ (falling back to the bundled data/example/ if none are found)
% and makes the whole-tissue (outside) view and/or the cross-section
% (ring) view, independently controlled by flag_tissue/
% flag_cross_section below -- run either one alone, or both. There is
% no true concurrency between them (this MATLAB install has no
% Parallel Computing Toolbox, so real parallel figure rendering isn't
% practical here): with both flags on they still run one after the
% other in this one script, but each is now independent -- turning
% tissue off means cross-section starts immediately, with no wait on a
% video you didn't ask for.
%
% Both go through the SAME video-making code path regardless of how
% many snapshots are selected: given an array of many timesteps it
% writes a real video; given an array of exactly one timestep it still
% "makes the video", which in that case is just that one plot, saved
% the same way.
%
% This script locates itself (rather than relying on the current
% working directory, which e.g. MATLAB's run() changes to the script's
% own folder) so it works regardless of how/from where it is invoked.

clear; clc; close all;

this_dir     = fileparts(mfilename('fullpath'));  % .../matlab
project_root = fileparts(this_dir);               % one level up
data_dir     = fullfile(project_root, 'data');
addpath(this_dir);

%% ---------------- user-configurable flags -----------------------------
flag_tissue        = true;           % make the whole-tissue (outside) view
flag_cross_section = true;          % make the cross-section (ring) view
                                      % (independent of flag_tissue -- run
                                      % either one alone, or both)

colorby = 'nsides';                  % 'nsides' | 'volume' | 'area' -- tissue colouring

cut_axis  = 'y';                     % cross-section cutting axis: 'x' | 'y' | 'z'
cut_value = 0.0;                     % cross-section cutting plane offset

% Which saved snapshots go into the video, given directly as an array
% of SAVED iteration numbers (the it_dumps cadence) -- not raw
% simulation steps, not a plain index into the file list. A scalar
% (a single iteration) works exactly like a 1-element array: the same
% code path still "makes the video", which in that case is just a
% single plot. Leave empty ([]) to use every available snapshot.
%
%   video_its = 1000:100:5000;   % that whole range
%   video_its = 5000;            % just one snapshot -> a single plot
%   video_its = [];               % every saved snapshot

video_its = 100000;

video_fps = 4;                     % frames per second for saved videos
video_dir = fullfile(project_root, 'videos');
%% ------------------------------------------------------------------------

[files, its] = list_snapshots(data_dir);
if isempty(files)
    % no simulation run yet (data/ is empty/gitignored) -- fall back to
    % the small curated example dataset shipped in the repo
    fprintf('No snapshots in %s -- using the bundled data/example/ instead.\n', data_dir);
    data_dir = fullfile(project_root, 'data', 'example');
    [files, its] = list_snapshots(data_dir);
end
if isempty(files)
    error('main_movie:nofiles', ...
          'No snapshots found in %s or its example/ subfolder -- run the simulation first.', ...
          fullfile(project_root, 'data'));
end
meta = read_mesh_meta(fullfile(data_dir, 'mesh_meta.txt'));

fprintf('Found %d snapshots (it = %d .. %d). R_apical=%.3g  R_basal=%.3g\n', ...
        numel(files), min(its), max(its), meta.R_apical, meta.R_basal);

[video_files, video_its_used] = select_snapshots(files, its, video_its);
fprintf('Selected %d frame(s): it = %s\n', numel(video_files), mat2str(video_its_used));

if ~exist(video_dir, 'dir')
    mkdir(video_dir);
end

if ~flag_tissue && ~flag_cross_section
    warning('main_movie:nothingtodo', 'Both flag_tissue and flag_cross_section are off -- nothing to do.');
end

if flag_tissue
    outfile = fullfile(video_dir, sprintf('tissue_%s.avi', colorby));
    make_tissue_video(video_files, colorby, outfile, video_fps);
end

if flag_cross_section
    outfile = fullfile(video_dir, sprintf('cross_section_%s_%s.avi', cut_axis, colorby));
    make_cross_section_video(video_files, cut_axis, cut_value, outfile, video_fps, colorby);
end
