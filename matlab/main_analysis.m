% main_analysis.m
% Example driver: reads the tissue snapshots written by the Fortran
% code into ../data/ and produces the two standard plots:
%   1. the whole tissue (standard 3D view)
%   2. a cutaway cross-section revealing the ring / hollow interior
%      alongside the cell shapes
%
% Run this script with the working directory set to the project root
% (the one containing data/), or edit data_dir below.

data_dir = 'data';
addpath(pwd);

meta  = read_mesh_meta(fullfile(data_dir, 'mesh_meta.txt'));
files = list_snapshots(data_dir);
if isempty(files)
    error('main_analysis:nofiles', 'No snapshots found in %s -- run the simulation first.', data_dir);
end

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
