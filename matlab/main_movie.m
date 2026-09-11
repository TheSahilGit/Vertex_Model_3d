% main_movie.m
% Driver: reads the tissue snapshots written by the Fortran code into
% data/ (falling back to the bundled data/example/ if none are found)
% and makes ONE combined video (make_combined_video.m) with three
% side-by-side views of every frame, all in one figure:
%   1. the whole tissue (outside view)
%   2. the LATITUDE ring through a chosen cell/location
%   3. the LONGITUDE ring through that same cell/location
% (see ring_planes_for_location.m for what "latitude"/"longitude" mean
% here: z is treated as the polar axis by convention).
%
% This script locates itself (rather than relying on the current
% working directory, which e.g. MATLAB's run() changes to the script's
% own folder) so it works regardless of how/from where it is invoked.

clear; clc; close all;

this_dir     = fileparts(mfilename('fullpath'));  % .../matlab
project_root = fileparts(this_dir);               % one level up
data_dir     = fullfile(project_root, 'data');
addpath(this_dir);

%% ---------------- user-configurable settings ---------------------------
colorby = 'nsides';                  % 'nsides' | 'volume' | 'area' -- shared
                                      % colouring for all three panels

ref_location = 1;                    % which cell/location the latitude and
                                      % longitude rings are cut through --
                                      % EITHER a cell index (e.g. 1) OR an
                                      % explicit [x y z] point (e.g.
                                      % [5 0 9]). Resolved once, from the
                                      % FIRST selected frame, into two
                                      % fixed cutting planes reused for
                                      % every later frame (see
                                      % ring_planes_for_location.m).

% Which saved snapshots go into the video, given directly as an array
% of SAVED iteration numbers (the it_dumps cadence) -- not raw
% simulation steps, not a plain index into the file list. A scalar
% (a single iteration) works exactly like a 1-element array: the same
% code path still "makes the video", which in that case is just a
% single (3-panel) plot. Leave empty ([]) to use every available snapshot.
%
%   video_its = 1000:100:5000;   % that whole range
%   video_its = 5000;            % just one snapshot -> a single plot
%   video_its = [];               % every saved snapshot

video_its = 100:100:10000;

video_fps = 4;                     % frames per second for the saved video
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

outfile = fullfile(video_dir, sprintf('combined_%s.avi', colorby));
make_combined_video(video_files, ref_location, colorby, outfile, video_fps);
