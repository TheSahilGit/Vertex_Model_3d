% main_analysis.m
% Driver: reads the tissue snapshots written by the Fortran code into
% data/ (falling back to the bundled data/example/ if none are found)
% and produces exactly two things, each controlled by the flags below:
%   - the whole-tissue (outside) view, always made
%   - the cutaway cross-section view, only if flag_cross_section is on
%
% Both go through the SAME video-making code path (make_tissue_video.m
% / make_cross_section_video.m) regardless of how many snapshots are
% selected: given an array of many timesteps it writes a real video;
% given an array of exactly one timestep it still "makes the video",
% which in that case is just that one plot, saved the same way. So
% there is no separate static-preview step duplicating figures --
% you get one tissue figure, and (if enabled) one cross-section
% figure, not one-per-colouring-mode-plus-a-preview.
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
colorby = 'nsides';                 % 'nsides' | 'volume' | 'area' -- tissue colouring

flag_cross_section = false;        % also make the cross-section video/plot
cut_axis  = 'y';                    % cross-section cutting axis: 'x' | 'y' | 'z'
cut_value = 0.0;                    % cross-section cutting plane offset

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

video_its = 1000:100:10000;

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
    error('main_analysis:nofiles', ...
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

% ---- the one tissue figure ----
outfile = fullfile(video_dir, sprintf('tissue_%s.avi', colorby));
make_tissue_video(video_files, colorby, outfile, video_fps);

% ---- the one, flagged, cross-section figure ----
if flag_cross_section
    outfile = fullfile(video_dir, sprintf('cross_section_%s.avi', cut_axis)); %#ok<UNRCH>
    make_cross_section_video(video_files, cut_axis, cut_value, outfile, video_fps);
end
