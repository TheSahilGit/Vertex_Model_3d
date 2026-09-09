% main_analysis.m
% Example driver: reads the tissue snapshots written by the Fortran
% code into ../data/ and produces the two standard plots:
%   1. the whole tissue (standard 3D view)
%   2. a cutaway cross-section revealing the ring / hollow interior
%      alongside the cell shapes
%
% This script locates itself (rather than relying on the current
% working directory, which e.g. MATLAB's run() changes to the script's
% own folder) so it works regardless of how/from where it is invoked.

this_dir     = fileparts(mfilename('fullpath'));  % .../matlab
project_root = fileparts(this_dir);               % one level up
data_dir     = fullfile(project_root, 'data');
addpath(this_dir);

files = list_snapshots(data_dir);
if isempty(files)
    % no simulation run yet (data/ is empty/gitignored) -- fall back to
    % the small curated example dataset shipped in the repo
    fprintf('No snapshots in %s -- using the bundled data/example/ instead.\n', data_dir);
    data_dir = fullfile(project_root, 'data', 'example');
    files = list_snapshots(data_dir);
end
if isempty(files)
    error('main_analysis:nofiles', ...
          'No snapshots found in %s or its example/ subfolder -- run the simulation first.', ...
          fullfile(project_root, 'data'));
end
meta = read_mesh_meta(fullfile(data_dir, 'mesh_meta.txt'));

fprintf('Found %d snapshots. R_apical=%.3g  R_basal=%.3g\n', ...
        numel(files), meta.R_apical, meta.R_basal);

% ---- last available snapshot: whole-tissue view + cross-section ----
S = read_vertex_snapshot(files{end});

plot_tissue_3d(S, 'nsides');
plot_cross_section(S, 'y', 0.0);

% ---- also show volume-strain and area-strain colouring, for a quick
% visual check that cells are near their target shape ----
plot_tissue_3d(S, 'volume');
plot_tissue_3d(S, 'area');

% ---- optional: step through every snapshot to watch the tissue
% evolve (uncomment to use) ----
% for k = 1:numel(files)
%     Sk = read_vertex_snapshot(files{k});
%     plot_tissue_3d(Sk, 'nsides');
%     drawnow;
%     pause(0.05);
% end
