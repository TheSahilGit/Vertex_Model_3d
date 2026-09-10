% main_analysis.m
% Driver: reads the scalar diagnostics time series written by the
% Fortran code (data/diag_<it>.dat, one small file per dump, alongside
% each snap_<it>.dat -- see the header comment of src/io_module.f90 for
% the binary layout) and plots them all against time:
%   - total energy of the system
%   - lumen volume (the hollow interior, enclosed by the basal surface)
%   - outer surface area (total apical area)
%   - max force (largest |force| over every alive vertex)
%   - number of (alive) cells
%   - cumulative T1 event count
%   - cumulative T2 event count
%
% For the tissue/cross-section VIDEOS (not this scalar time series),
% see main_movie.m instead.
%
% This script locates itself (rather than relying on the current
% working directory, which e.g. MATLAB's run() changes to the script's
% own folder) so it works regardless of how/from where it is invoked.

clear; clc; close all;

this_dir     = fileparts(mfilename('fullpath'));  % .../matlab
project_root = fileparts(this_dir);               % one level up
data_dir     = fullfile(project_root, 'data');
addpath(this_dir);

FONT_SIZE = 34;

[files, its] = list_diagnostics(data_dir);
if isempty(files)
    % no simulation run yet (data/ is empty/gitignored) -- fall back to
    % the small curated example dataset shipped in the repo
    fprintf('No diagnostics in %s -- using the bundled data/example/ instead.\n', data_dir);
    data_dir = fullfile(project_root, 'data', 'example');
    [files, its] = list_diagnostics(data_dir);
end
if isempty(files)
    error('main_analysis:nofiles', ...
          'No diagnostics found in %s or its example/ subfolder -- run the simulation first.', ...
          fullfile(project_root, 'data'));
end

[its, order] = sort(its);
files = files(order);

n = numel(files);
time          = zeros(n, 1);
energy        = zeros(n, 1);
lumen_volume  = zeros(n, 1);
outer_area    = zeros(n, 1);
max_force     = zeros(n, 1);
n_cells       = zeros(n, 1);
cumulative_T1 = zeros(n, 1);
cumulative_T2 = zeros(n, 1);

for k = 1:n
    D = read_diagnostics(files{k});
    time(k)          = D.time;
    energy(k)        = D.energy;
    lumen_volume(k)  = D.lumen_volume;
    outer_area(k)    = D.outer_area;
    max_force(k)     = D.max_force;
    n_cells(k)       = D.n_cells;
    cumulative_T1(k) = D.cumulative_T1;
    cumulative_T2(k) = D.cumulative_T2;
end

fprintf('Found %d diagnostics dumps (it = %d .. %d)\n', n, its(1), its(end));

series  = {energy, lumen_volume, outer_area, max_force, n_cells, cumulative_T1, cumulative_T2};
titles  = {'Total energy', 'Lumen volume (hollow interior)', 'Outer surface area', ...
           'Max |force| over all vertices', 'Number of cells', ...
           'Cumulative T1 events', 'Cumulative T2 events'};
ylabels = {'E', 'V_{lumen}', 'A_{outer}', 'F_{max}', 'N_{cell}', 'N_{T1}', 'N_{T2}'};

fig = figure('Color', 'w', 'Position', [100 100 1500 950]);
for idx = 1:numel(series)
    subplot(3, 3, idx);
    plot(time, series{idx}, '-', 'LineWidth', 6, 'Color', [0.15 0.35 0.75]);
    xlabel('t', 'FontSize', FONT_SIZE);
    ylabel(ylabels{idx}, 'FontSize', FONT_SIZE);
    %title(titles{idx}, 'FontSize', FONT_SIZE);
    set(gca, 'FontSize', FONT_SIZE , 'LineWidth', 4, 'FontName', 'sans');
    grid off
end
%sgtitle('Simulation diagnostics vs time', 'FontSize', FONT_SIZE + 4, 'FontWeight', 'bold');
