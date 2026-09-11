% main_movie.m
% Driver: reads the tissue snapshots written by the Fortran code into
% data/ (falling back to the bundled data/example/ if none are found)
% and makes up to THREE independent videos, each its own figure/window:
%   1. the whole tissue (outside view)                 -- flag_tissue
%   2. the LATITUDE ring through a chosen cell/location -- flag_latitude
%   3. the LONGITUDE ring through that same location    -- flag_longitude
% (see ring_planes_for_location.m for what "latitude"/"longitude" mean
% here: z is treated as the polar axis by convention). The three flags
% are independent -- run any subset of them, in any combination; there
% is no shared figure between them, so rotating one interactively can
% never affect another.
%
% Each goes through the SAME video-making code path regardless of how
% many snapshots are selected: given an array of many timesteps it
% writes a real video; given an array of exactly one timestep it still
% "makes the video", which in that case is just that one plot, saved
% the same way.
%
% With more than one flag on, they run one after another, not at the
% same time -- checked directly: this MATLAB install is licensed for
% Parallel Computing Toolbox but does not actually have it installed
% (parpool/gcp are undefined here), so there is no parfor/parfeval
% available to run them concurrently within one MATLAB session.
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
flag_tissue     = true;    % make the whole-tissue (outside) view
flag_latitude   = true;    % make the latitude-ring cross-section view
flag_longitude  = true;    % make the longitude-ring cross-section view
                            % (all three independent -- run any subset)

colorby = 'shapefactor';         % colouring for all three (individual cells
                             % shown as one solid colour each on the
                             % rings, too) -- see tissue_color_values.m
                             % for the full list and exact definitions:
                             %   'nsides'                    number of sides
                             %   'volume' / 'volume_abs'     V/V0-1 / V
                             %   'area' / 'apical_area'      Aapi/A0-1 / Aapi
                             %   'basal_area' / 'basal_area_abs'   Abas/A0bas-1 / Abas
                             %   'lateral_area'              cell-cell wall area
                             %   'total_area'                apical+basal+lateral
                             %   'shapefactor'               total_area / V^(2/3)
colormapName = 'parula';    % any built-in MATLAB colormap name, e.g.
                             % 'turbo', 'jet', 'hot', 'cool', 'copper', 'bone'

ref_location = 10;            % which cell/location the latitude and
                              % longitude rings are cut through -- EITHER
                              % a cell index (e.g. 1) OR an explicit
                              % [x y z] point (e.g. [5 0 9]). Resolved
                              % once, from the FIRST selected frame, into
                              % two fixed cutting planes reused for every
                              % later frame (see ring_planes_for_location.m).

latitude_view  = [-124.6069609926535,  33.152764627008267];  % [az el],
longitude_view = [  -9.195907331898637, 10.77613463554659];  % found by hand
                              % (rotate the figure, then read back with
                              % [az,el]=view(gca)) -- used exactly as
                              % given instead of the automatic tilt.

% Which saved snapshots go into the video(s), given directly as an
% array of SAVED iteration numbers (the it_dumps cadence) -- not raw
% simulation steps, not a plain index into the file list. A scalar
% (a single iteration) works exactly like a 1-element array: the same
% code path still "makes the video", which in that case is just a
% single plot. Leave empty ([]) to use every available snapshot.
%
%   video_its = 1000:100:5000;   % that whole range
%   video_its = 5000;            % just one snapshot -> a single plot
%   video_its = [];               % every saved snapshot

video_its =100:1000:100000;

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

if ~flag_tissue && ~flag_latitude && ~flag_longitude
    warning('main_movie:nothingtodo', 'All three flags are off -- nothing to do.');
end

if flag_tissue
    outfile = fullfile(video_dir, sprintf('tissue_%s.avi', colorby));
    make_tissue_video(video_files, colorby, outfile, video_fps, colormapName);
end

if flag_latitude || flag_longitude
    S1 = read_vertex_snapshot(video_files{1});
    [latPlane, lonPlane, refPoint] = ring_planes_for_location(S1, ref_location);
    fprintf('Ring reference point = [%.4g %.4g %.4g]\n', refPoint);

    if flag_latitude
        outfile = fullfile(video_dir, sprintf('latitude_%s.avi', colorby));
        make_cross_section_video(video_files, latPlane, outfile, video_fps, colorby, ...
                                  'Latitude ring', latitude_view, colormapName);
    end
    if flag_longitude
        outfile = fullfile(video_dir, sprintf('longitude_%s.avi', colorby));
        make_cross_section_video(video_files, lonPlane, outfile, video_fps, colorby, ...
                                  'Longitude ring', longitude_view, colormapName);
    end
end
