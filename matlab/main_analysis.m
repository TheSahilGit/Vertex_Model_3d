% main_analysis.m
% Driver: reads the tissue snapshots written by the Fortran code into
% data/ (falling back to the bundled data/example/ if none are found)
% and, controlled by the flags below, produces:
%   - a quick static preview of the latest snapshot: whole-tissue view
%     (per enabled colouring) + cutaway cross-section
%   - a video stepping through a chosen range of snapshots, for the
%     whole-tissue view (one video per enabled colouring)
%
% This script locates itself (rather than relying on the current
% working directory, which e.g. MATLAB's run() changes to the script's
% own folder) so it works regardless of how/from where it is invoked.

this_dir     = fileparts(mfilename('fullpath'));  % .../matlab
project_root = fileparts(this_dir);               % one level up
data_dir     = fullfile(project_root, 'data');
addpath(this_dir);

%% ---------------- user-configurable flags -----------------------------
flag_video_tissue        = true;    % save a video: whole-tissue (outside) view
flag_video_cross_section = false;   % cross-section video: revisit later

flag_color_nsides = true;          % include "number of sides" colouring (tissue view/video)
flag_color_volume = true;          % include "volume strain" colouring (tissue view/video)
flag_color_area   = true;          % include "apical area strain" colouring (tissue view/video)

cut_axis  = 'y';                   % cross-section cutting axis: 'x' | 'y' | 'z'
cut_value = 0.0;                   % cross-section cutting plane offset

% Which saved snapshots go into the tissue video, expressed in terms
% of the SAVED iteration numbers (the it_dumps cadence), not raw
% simulation steps and not a plain index into the file list:
video_it_start  = -inf;   % first saved iteration to include (-inf = from the very first)
video_it_end    = inf;    % last saved iteration to include  ( inf = up to the very last)
video_it_stride = 1;      % use every Nth saved snapshot in that range (1 = all of them)

video_fps = 8;                     % frames per second for saved videos
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

colorby_list = {};
if flag_color_nsides, colorby_list{end+1} = 'nsides'; end %#ok<UNRCH>
if flag_color_volume, colorby_list{end+1} = 'volume'; end
if flag_color_area,   colorby_list{end+1} = 'area';   end

% ---- quick static preview of the latest snapshot ----
S = read_vertex_snapshot(files{end});
for k = 1:numel(colorby_list)
    plot_tissue_3d(S, colorby_list{k});
end
plot_cross_section(S, cut_axis, cut_value);

% ---- tissue video across the selected snapshot range ----
if flag_video_tissue
    [video_files, video_its] = select_snapshots(files, its, video_it_start, video_it_end, video_it_stride);
    fprintf('Tissue video: %d frames selected (it = %d .. %d, stride %d)\n', ...
            numel(video_files), min(video_its), max(video_its), video_it_stride);

    if ~exist(video_dir, 'dir')
        mkdir(video_dir);
    end
    for k = 1:numel(colorby_list)
        cb = colorby_list{k};
        outfile = fullfile(video_dir, sprintf('tissue_%s.avi', cb));
        make_tissue_video(video_files, cb, outfile, video_fps);
    end
end

% ---- cross-section video: revisit later ----
if flag_video_cross_section
    [video_files, video_its] = select_snapshots(files, its, video_it_start, video_it_end, video_it_stride); %#ok<UNRCH>
    if ~exist(video_dir, 'dir')
        mkdir(video_dir);
    end
    outfile = fullfile(video_dir, sprintf('cross_section_%s.avi', cut_axis));
    make_cross_section_video(video_files, cut_axis, cut_value, outfile, video_fps);
end
