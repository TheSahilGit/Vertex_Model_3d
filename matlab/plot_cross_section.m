function fig = plot_cross_section(S, plane, varargin)
% PLOT_CROSS_SECTION  A single RING of cells -- the ones whose apical
% polygon straddles a cutting plane -- rendered as complete polyhedra
% (apical + lateral + basal faces, all edges visible), each cell
% coloured as one solid colour by a per-cell metric (like
% plot_tissue_3d.m) so individual cells are easy to tell apart, and
% viewed face-on to the cutting plane. Even though it's a genuine 3D
% structure, viewed this way it reads almost as a flat 2D strip of
% polygons: exactly the ring of cells a physical cut through the
% tissue would reveal, without the foreshortening/distortion of
% rendering the whole 3D shell from an oblique angle.
%
% PLANE identifies the cutting plane, in either of two forms:
%
%   plot_cross_section(S, 'y', 0)          % legacy form: cutaxis/cutvalue,
%                                           % 'x'|'y'|'z' plus a scalar --
%                                           % the plane {axis = value}
%   plot_cross_section(S, planeStruct)     % general form: a struct with
%                                           % .normal (3x1) and .point (3x1),
%                                           % the plane {r : dot(r-point,normal)=0}
%                                           % -- e.g. from ring_planes_for_location.m,
%                                           % for a plane at an arbitrary angle
%                                           % (not just axis-aligned)
%
%   plot_cross_section(S)                        % cut at y = 0, colour by volume
%   plot_cross_section(S, 'z', 2.0)              % cut at z = 2.0
%   plot_cross_section(S, latPlane, 'ColorBy', 'nsides')
%   plot_cross_section(S, lonPlane, 'Figure', fig, 'CLim', [cmin cmax])
%       draw into an existing figure (cleared first), with a fixed
%       colour range instead of auto-scaling to this frame's data --
%       used by make_cross_section_video.m for a constant window/frame
%       size and a constant colour scale across a whole video.
%   plot_cross_section(S, latPlane, 'Axes', ax, 'CLim', [cmin cmax])
%       draw into an existing AXES (e.g. one subplot of a bigger
%       figure) instead of owning a whole figure -- used by
%       make_combined_video.m. Takes precedence over 'Figure'.
%
% S is a struct from read_vertex_snapshot.m.

if nargin < 2 || isempty(plane), plane = 'y'; end

% ---- resolve PLANE into a unit normal + a point on the plane ----
if ischar(plane) || isstring(plane)
    cutaxis = char(plane);
    if nargin < 3 || isempty(varargin) || ~isnumeric(varargin{1})
        cutvalue = 0.0;
        rest = varargin;
    else
        cutvalue = varargin{1};
        rest = varargin(2:end);
    end
    switch lower(cutaxis)
        case 'x', planeNormal = [1;0;0];
        case 'y', planeNormal = [0;1;0];
        case 'z', planeNormal = [0;0;1];
        otherwise
            error('plot_cross_section:cutaxis', 'cutaxis must be ''x'', ''y'' or ''z''');
    end
    planePoint = planeNormal * cutvalue;
else
    planeNormal = plane.normal(:) / norm(plane.normal);
    planePoint  = plane.point(:);
    rest = varargin;
end

p = inputParser;
p.addParameter('Figure', []);
p.addParameter('Axes', []);
p.addParameter('ColorBy', 'volume');
p.addParameter('CLim', []);
p.parse(rest{:});
fig_in   = p.Results.Figure;
ax_in    = p.Results.Axes;
colorby  = p.Results.ColorBy;
clim_in  = p.Results.CLim;

FONT_SIZE  = 26;
TITLE_SIZE = 30;

% ---- per-cell colour value for every alive cell (same definition as
% plot_tissue_3d.m), looked up per ring cell below ----
[cval_alive, cblabel] = tissue_color_values(S, colorby);
cval_map = zeros(S.n_cell, 1);
cval_map(S.cell_alive) = cval_alive;

% ---- the ring: alive cells whose apical polygon straddles the plane ----
idx_alive = find(S.cell_alive);
ring_idx = zeros(numel(idx_alive), 1);
nring = 0;
for ii = 1:numel(idx_alive)
    ic = idx_alive(ii);
    n = S.cell_n(ic);
    vids = S.cell_vlist(ic, 1:n);
    d = (S.r_api(vids, :) - planePoint') * planeNormal;
    if any(d > 0) && any(d < 0)
        nring = nring + 1;
        ring_idx(nring) = ic;
    end
end
ring_idx = ring_idx(1:nring);

if nring == 0
    warning('plot_cross_section:noring', ...
        'No cells straddle the requested plane; nothing to plot.');
end

% ---- explicit triangles for every face (apical/basal fans from each
% cell's own centroid, lateral quads split into 2 triangles), rather
% than handing MATLAB N-gon/quad faces to triangulate internally via
% patch's own Delaunay-based tessellation. That was observed to
% occasionally go numerically unstable ("Enforcing Delaunay constraints
% leads to unstable repeated constraint intersections") for a
% particular cell at a particular frame of an otherwise fine long
% video -- this ring is small enough (a few dozen cells) that doing
% the triangulation ourselves, exactly like Force.f90 already does for
% the real physics, is cheap and removes the dependency entirely.
%
% Every triangle from a given cell carries that cell's colour value as
% its FaceVertexCData entry (one row per face), so with FaceColor
% 'flat' each cell -- apical cap, basal cap and all its lateral walls
% alike -- renders as one solid, distinct colour.
%
% Face edges are drawn separately as plain line segments (not via each
% patch's own EdgeColor), since the fan/split triangulation would
% otherwise show its internal spoke/diagonal lines as spurious extra
% edges; only each face's true polygon boundary is drawn.
nvert = size(S.r_api, 1);   % basal vertex k sits at row nvert+k in combinedV
% combinedV row layout: [r_api (1..nvert); r_bas (1..nvert);
%                        apical centroids (1..nring); basal centroids (1..nring)]
apiC_base = 2*nvert;
basC_base = 2*nvert + nring;

Tapi = zeros(0, 3); Capi = zeros(0, 1);
Tbas = zeros(0, 3); Cbas = zeros(0, 1);
Tlat = zeros(0, 3); Clat = zeros(0, 1);
apiCentroids = zeros(nring, 3);
basCentroids = zeros(nring, 3);
edge_x = []; edge_y = []; edge_z = [];

for k = 1:nring
    ic = ring_idx(k);
    n = S.cell_n(ic);
    vids = S.cell_vlist(ic, 1:n);
    api_ids = vids;
    bas_ids = vids + nvert;
    ca_id = apiC_base + k;
    cb_id = basC_base + k;
    cval = cval_map(ic);

    apiCentroids(k, :) = mean(S.r_api(vids, :), 1);
    basCentroids(k, :) = mean(S.r_bas(vids, :), 1);

    for e = 1:n
        e2 = mod(e, n) + 1;
        Tapi(end+1, :) = [ca_id, api_ids(e), api_ids(e2)]; Capi(end+1, 1) = cval; %#ok<AGROW>
        Tbas(end+1, :) = [cb_id, bas_ids(e2), bas_ids(e)]; Cbas(end+1, 1) = cval; %#ok<AGROW> % reversed: outward-facing basal cap

        a1 = api_ids(e); a2 = api_ids(e2);
        b1 = bas_ids(e); b2 = bas_ids(e2);
        Tlat(end+1, :) = [a1, a2, b2]; Clat(end+1, 1) = cval; %#ok<AGROW>
        Tlat(end+1, :) = [a1, b2, b1]; Clat(end+1, 1) = cval; %#ok<AGROW>

        edge_x = [edge_x; S.r_api(a1,1); S.r_api(a2,1); NaN; ...   % apical edge
                           S.r_bas(a1,1); S.r_bas(a2,1); NaN; ...   % basal edge
                           S.r_api(a1,1); S.r_bas(a1,1); NaN];      % lateral (vertical) edge
        edge_y = [edge_y; S.r_api(a1,2); S.r_api(a2,2); NaN; ...
                           S.r_bas(a1,2); S.r_bas(a2,2); NaN; ...
                           S.r_api(a1,2); S.r_bas(a1,2); NaN];
        edge_z = [edge_z; S.r_api(a1,3); S.r_api(a2,3); NaN; ...
                           S.r_bas(a1,3); S.r_bas(a2,3); NaN; ...
                           S.r_api(a1,3); S.r_bas(a1,3); NaN];
    end
end

combinedV = [S.r_api; S.r_bas; apiCentroids; basCentroids];

if ~isempty(ax_in)
    ax = ax_in;
    cla(ax);
    fig = ancestor(ax, 'figure');
else
    if isempty(fig_in)
        fig = figure('Color', 'w');
    else
        fig = fig_in;
        figure(fig);
        clf(fig);
    end
    ax = axes(fig);
end

patch(ax, 'Faces', Tapi, 'Vertices', combinedV, 'FaceVertexCData', Capi, 'FaceColor', 'flat', 'EdgeColor', 'none');
hold(ax, 'on')
patch(ax, 'Faces', Tbas, 'Vertices', combinedV, 'FaceVertexCData', Cbas, 'FaceColor', 'flat', 'EdgeColor', 'none');
patch(ax, 'Faces', Tlat, 'Vertices', combinedV, 'FaceVertexCData', Clat, 'FaceColor', 'flat', 'EdgeColor', 'none');
plot3(ax, edge_x, edge_y, edge_z, 'Color', [0.05 0.05 0.05], 'LineWidth', 1.2);

axis(ax, 'equal'); axis(ax, 'off');

colormap(ax, parula)
if isempty(clim_in)
    if max(cval_alive) > min(cval_alive)
        clim(ax, [min(cval_alive), max(cval_alive)]);
    end
else
    clim(ax, clim_in);
end
cb = colorbar(ax);
cb.Label.String = cblabel;
cb.FontSize = FONT_SIZE;

% Fixed reference frame (constant R_apical, not this frame's own data
% extent) -- see plot_tissue_3d.m for why: it prevents MATLAB's
% auto-zoom from re-fitting to the tissue's own (genuinely, slowly
% shrinking) size, which otherwise showed up as an occasional sudden
% jump in rendered size partway through an otherwise perfectly smooth
% video. (This box just keeps the DATA limits fixed; the actual
% on-screen framing below is set explicitly and does not depend on it.)
R = S.R_apical * 1.05;
xlim(ax, [-R R]); ylim(ax, [-R R]); zlim(ax, [-R R]);

% Face-on to the cutting plane, with an on-screen framing that stays
% CONSISTENT regardless of which direction planeNormal happens to
% point. view(ax, v) alone would auto-fit the camera distance/angle to
% whatever the fixed [-R,R]^3 axis box's SILHOUETTE looks like from
% that specific direction -- and that silhouette is wider when viewed
% along an oblique in-plane direction (e.g. a longitude ring at an
% arbitrary azimuthal angle, up to R*sqrt(2) across) than dead-on
% along a coordinate axis (e.g. the latitude ring, looking straight
% down z, exactly R across) -- so two rings of the SAME physical
% radius R would render at visibly different on-screen sizes purely
% because of their different cut orientation. Setting camera position/
% target/view-angle explicitly (a fixed distance and view angle
% computed from R alone, independent of direction) removes that
% dependence: every ring, whatever its cutting plane, fills the same
% fraction of its panel.
camdist = 4 * S.R_apical;
camtarget(ax, [0 0 0]);
campos(ax, planeNormal' * camdist);
upvec = [0 0 1];
if abs(dot(upvec, planeNormal)) > 0.9
    upvec = [0 1 0];   % avoid a degenerate/near-parallel up vector when
                        % looking nearly straight down the polar axis
end
camup(ax, upvec);
camva(ax, 2 * atand(R / camdist));

camlight(ax, 'headlight'); lighting(ax, 'gouraud'); material(ax, 'dull')
ax.CameraPositionMode  = 'manual';
ax.CameraTargetMode    = 'manual';
ax.CameraUpVectorMode  = 'manual';
ax.CameraViewAngleMode = 'manual';  % freeze zoom AFTER Position/limits are final
ax.FontSize = FONT_SIZE;
ax.Toolbar.Visible = 'on';  % 'off' to avoid it showing up in exported/captured frames

title(ax, sprintf('t = %.4g', S.time), 'FontSize', TITLE_SIZE);
end
