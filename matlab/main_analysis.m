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
%   - cumulative T4 (crowding-induced extrusion) event count, and its
%     apical-ward / basal-ward / ambiguous breakdown
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
cumulative_T4           = zeros(n, 1);
cumulative_T4_apical    = zeros(n, 1);
cumulative_T4_basal     = zeros(n, 1);
cumulative_T4_ambiguous = zeros(n, 1);

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
    cumulative_T4(k)           = D.cumulative_T4;
    cumulative_T4_apical(k)    = D.cumulative_T4_apical;
    cumulative_T4_basal(k)     = D.cumulative_T4_basal;
    cumulative_T4_ambiguous(k) = D.cumulative_T4_ambiguous;
end

fprintf('Found %d diagnostics dumps (it = %d .. %d)\n', n, its(1), its(end));

series  = {energy, lumen_volume, outer_area, max_force, n_cells, cumulative_T1, cumulative_T2, cumulative_T4};
titles  = {'Total energy', 'Lumen volume (hollow interior)', 'Outer surface area', ...
           'Max |force| over all vertices', 'Number of cells', ...
           'Cumulative T1 events', 'Cumulative T2 events', 'Cumulative T4 events'};
ylabels = {'E', 'V_{lumen}', 'A_{outer}', 'F_{max}', 'N_{cell}', 'N_{T1}', 'N_{T2}', 'N_{T4}'};

fig = figure('Color', 'w', 'Position', [100 100 1500 950]);
for idx = 1:numel(series)
    subplot(3, 3, idx);
    % if strcmp(ylabels{idx}, 'E') || strcmp(ylabels{idx}, 'F_{max}') | strcmp(ylabels{idx}, 'V_{lumen}')
    %     loglog(time, series{idx}, '-', 'LineWidth', 6, 'Color', [0.15 0.35 0.75]);
    % else
    % plot(time, series{idx}, '-', 'LineWidth', 6, 'Color', [0.15 0.35 0.75]);
    % end
    loglog(time, series{idx}, '-', 'LineWidth', 6, 'Color', [0.15 0.35 0.75]);

    xlabel('t', 'FontSize', FONT_SIZE);
    ylabel(ylabels{idx}, 'FontSize', FONT_SIZE);
    %title(titles{idx}, 'FontSize', FONT_SIZE);
    set(gca, 'FontSize', FONT_SIZE , 'LineWidth', 4, 'FontName', 'sans');
    grid off
end

% Slot 9: T4's apical-ward/basal-ward/ambiguous breakdown -- three lines,
% not a single series, so handled separately from the generic loop above.
subplot(3, 3, 9);
loglog(time, cumulative_T4_apical, '-', 'LineWidth', 6, 'Color', [0.85 0.30 0.30]);
hold on
loglog(time, cumulative_T4_basal, '-', 'LineWidth', 6, 'Color', [0.30 0.55 0.85]);
loglog(time, cumulative_T4_ambiguous, '-', 'LineWidth', 6, 'Color', [0.55 0.55 0.55]);
hold off
xlabel('t', 'FontSize', FONT_SIZE);
ylabel('N_{T4}', 'FontSize', FONT_SIZE);
legend({'apical', 'basal', 'ambiguous'}, 'FontSize', FONT_SIZE*0.5, 'Location', 'best');
set(gca, 'FontSize', FONT_SIZE, 'LineWidth', 4, 'FontName', 'sans');
grid off
%sgtitle('Simulation diagnostics vs time', 'FontSize', FONT_SIZE + 4, 'FontWeight', 'bold');
