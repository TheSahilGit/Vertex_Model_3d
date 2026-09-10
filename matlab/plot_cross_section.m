function fig = plot_cross_section(S, cutaxis, cutvalue, varargin)
% PLOT_CROSS_SECTION  Cutaway view showing the hollow interior, the
% ring (shell) cross-section, and the individual cell shapes.
%
% A cutting plane {cutaxis = cutvalue} through (or near) the sphere
% centre is used to keep only the cells on one side of it (a "cutaway
% hemisphere"), so their apical AND basal faces are both rendered,
% making the shell thickness and the hollow interior directly visible.
% A reference disk (radius R_basal, the exact intersection of the
% basal shell with the cutting plane) is drawn flat in that plane to
% unambiguously mark the hollow interior even where no cell happens to
% sit exactly on the plane.
%
%   plot_cross_section(S)                  % cut at y = 0
%   plot_cross_section(S, 'z', 2.0)        % cut at z = 2.0
%   plot_cross_section(S, 'y', 0, 'Figure', fig)
%       draw into an existing figure (cleared first) instead of
%       creating a new one -- used by make_cross_section_video.m for a
%       constant window/frame size across a whole video.
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2 || isempty(cutaxis),  cutaxis  = 'y'; end
if nargin < 3 || isempty(cutvalue), cutvalue = 0.0; end
p = inputParser;
p.addParameter('Figure', []);
p.parse(varargin{:});
fig_in = p.Results.Figure;

FONT_SIZE  = 26;  % ~2.5x the MATLAB default (10), per user request
TITLE_SIZE = 28;

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

if isempty(fig_in)
    fig = figure('Color', 'w');
else
    fig = fig_in;
    figure(fig);
    clf(fig);
end
hold on

% ---- flat reference disk at the basal radius, coloured white, in the
% cutting plane: from this oblique 3/4 viewpoint it sits just behind
% the real (curved) basal cell patches and "punches out" a clean blank
% hollow interior wherever no cell patch already covers it -- making
% the true hollow (r < R_basal) unambiguous even between cells. (An
% equivalent outer disk at R_apical was tried too, but at this oblique
% angle it is a flat disk sitting slightly proud of the curved apical
% shell and shows through the gaps between cells as odd flat shards;
% the real apical cell patches already mark the outer boundary fine
% without it.)
theta = linspace(0, 2*pi, 200)';
draw_disk(S.R_basal, cutvalue, axcol, [1 1 1], theta);

% ---- the actual cutaway cell shapes: apical (outer) and basal (inner)
% faces of every retained cell, so the shell thickness is explicit ----
h_api = patch('Faces', Fk, 'Vertices', S.r_api, 'FaceColor', [0.55 0.70 0.95], ...
      'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.5, 'FaceAlpha', 1.0);
h_bas = patch('Faces', Fk, 'Vertices', S.r_bas, 'FaceColor', [0.95 0.65 0.45], ...
      'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.5, 'FaceAlpha', 1.0);

% No 'vis3d' (which freezes the camera view angle for interactive
% rotation, not needed for a scripted, fixed-view plot/video): it was
% found to fight with an explicit axes-Position resize in
% plot_tissue_3d.m, unpredictably breaking the rendered size/title
% frame to frame -- see that file's comment for the full story.
axis equal off
% Flat, unlit colouring on purpose: with directional lighting, cells
% whose polygon is nearly edge-on to this viewing angle (e.g. ones
% that straddle the cut boundary, since cells are kept/dropped by
% apical-centroid side rather than true polygon clipping) get shaded
% almost black at grazing incidence, which reads as a rendering
% glitch. Flat faces keep the schematic honest and legible instead.

% Camera on the side the cut face's outward normal points to (the
% "kept" material is on the OTHER side, so this looks squarely into
% the cut), with a bit of extra azimuth/elevation for a 3/4, not
% perfectly flat-on, view. In MATLAB's view(az,el), az=0/90/180/270
% correspond to the camera sitting on the -y/+x/+y/-x axis
% respectively; since only the apical/basal CAPS are drawn here (no
% explicit lateral wall faces), a face-on view (az aligned with the
% cut axis, el~0) shows nothing -- verified empirically against real
% renders before picking these angles.
switch axcol
    case 1, view(110, 18);   % x-cut:  kept x<=cutvalue, cut normal +x
    case 2, view(200, 18);   % y-cut:  kept y<=cutvalue, cut normal +y
    case 3, view(35, 75);    % z-cut:  kept z<=cutvalue, cut normal +z (viewed from above)
end

set(gca, 'FontSize', FONT_SIZE);
% title(sprintf('Cross-section (%s = %.3g) at step %d  (t = %.4g)', ...
%       cutaxis, cutvalue, S.it, S.time), 'FontSize', TITLE_SIZE); %#ok<UNRCH>
legend([h_api, h_bas], {'apical (outer) faces', 'basal (inner) faces'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal', ...
       'FontSize', FONT_SIZE);
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
