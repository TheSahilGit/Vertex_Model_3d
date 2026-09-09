function fig = plot_cross_section(S, cutaxis, cutvalue)
% PLOT_CROSS_SECTION  Cutaway view showing the hollow interior, the
% ring (shell) cross-section, and the individual cell shapes.
%
% A cutting plane {cutaxis = cutvalue} through (or near) the sphere
% centre is used to keep only the cells on one side of it (a "cutaway
% hemisphere"), so their apical AND basal faces are both rendered,
% making the shell thickness and the hollow interior directly visible.
% Two reference circles (radius R_apical and R_basal, the exact
% intersection of the two shells with the cutting plane) are drawn
% flat in that plane, filled as an annulus, to unambiguously mark the
% ring and the hollow interior even where no cell happens to sit
% exactly on the plane.
%
%   plot_cross_section(S)                  % cut at y = 0
%   plot_cross_section(S, 'z', 2.0)        % cut at z = 2.0
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2 || isempty(cutaxis),  cutaxis  = 'y'; end
if nargin < 3 || isempty(cutvalue), cutvalue = 0.0; end

switch lower(cutaxis)
    case 'x', axcol = 1;
    case 'y', axcol = 2;
    case 'z', axcol = 3;
    otherwise
        error('plot_cross_section:cutaxis', 'cutaxis must be ''x'', ''y'' or ''z''');
end

[F, idx] = cell_faces_matrix(S);

% keep cells whose apical-face centroid lies on the "kept" side
keep = false(size(F, 1), 1);
for k = 1:size(F, 1)
    vids = F(k, ~isnan(F(k, :)));
    c = mean(S.r_api(vids, :), 1);
    keep(k) = c(axcol) <= cutvalue;
end
Fk = F(keep, :);

fig = figure('Color', 'w');
hold on

% ---- reference annulus in the cutting plane: outer disk then an
% inner disk drawn in the figure background colour "punches out" the
% hollow interior, showing the ring cross-section unambiguously ----
theta = linspace(0, 2*pi, 200)';
draw_disk(S.R_apical, cutvalue, axcol, [0.75 0.78 0.82], theta);
draw_disk(S.R_basal,  cutvalue, axcol, [1 1 1],          theta);

% ---- the actual cutaway cell shapes: apical (outer) and basal (inner)
% faces of every retained cell, so the shell thickness is explicit ----
h_api = patch('Faces', Fk, 'Vertices', S.r_api, 'FaceColor', [0.55 0.70 0.95], ...
      'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.5, 'FaceAlpha', 1.0);
h_bas = patch('Faces', Fk, 'Vertices', S.r_bas, 'FaceColor', [0.95 0.65 0.45], ...
      'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.5, 'FaceAlpha', 1.0);

axis equal vis3d off
camlight('headlight'); lighting gouraud; material dull

% look roughly along the cut normal so the cut face is fully visible
switch axcol
    case 1, view(-90, 10);
    case 2, view(0, 10);
    case 3, view(35, 75);
end

title(sprintf('Cross-section (%s = %.3g) at step %d  (t = %.4g)', ...
      cutaxis, cutvalue, S.it, S.time));
legend([h_api, h_bas], {'apical (outer) faces', 'basal (inner) faces'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal');
end

%==========================================================================
function draw_disk(R, cutvalue, axcol, faceColor, theta)
rho2 = R^2 - cutvalue^2;
if rho2 <= 0
    return  % cutting plane misses this shell entirely
end
rho = sqrt(rho2);
n = numel(theta);
P = zeros(n, 3);
switch axcol
    case 1
        P(:,1) = cutvalue; P(:,2) = rho*cos(theta); P(:,3) = rho*sin(theta);
    case 2
        P(:,2) = cutvalue; P(:,1) = rho*cos(theta); P(:,3) = rho*sin(theta);
    case 3
        P(:,3) = cutvalue; P(:,1) = rho*cos(theta); P(:,2) = rho*sin(theta);
end
patch('Faces', 1:n, 'Vertices', P, 'FaceColor', faceColor, ...
      'EdgeColor', 'none', 'FaceAlpha', 0.9);
end
