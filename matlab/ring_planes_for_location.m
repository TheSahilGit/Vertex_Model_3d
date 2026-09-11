function [latPlane, lonPlane, refPoint] = ring_planes_for_location(S, ref)
% RING_PLANES_FOR_LOCATION  Given a snapshot S and a reference location
% REF, return the two cutting planes through that location that
% plot_cross_section.m expects -- one giving the LATITUDE ring (the
% ring of cells at the same "height" along the polar axis) and one
% giving the LONGITUDE ring (the ring of cells at the same azimuthal
% angle -- the great-circle meridian plane through the poles) -- plus
% the resolved reference point itself.
%
% The polar axis is fixed as the z-axis by convention: this mesh has
% no intrinsically marked pole, and z is simply the axis every other
% plotting routine here already treats as "up" (plot_tissue_3d.m's
% default view, etc).
%
%   ref = 7                  % a cell index: uses that cell's OWN
%                             % apical-face centroid as the location
%   ref = [x y z]             % an explicit point in space instead
%
% Both REF forms are resolved ONCE (typically from a video's first
% frame -- see main_movie.m) into a single fixed reference
% point; the two planes returned are meant to be reused unchanged for
% every frame of a video, exactly like cut_axis/cut_value were fixed
% for a whole video before this. If the tissue moves or a tracked
% cell's own shape drifts, the rings drawn from these fixed planes
% simply follow whatever the tissue looks like there now -- they are
% NOT re-centred frame to frame.
%
% Each returned plane is a struct with fields .normal (3x1 unit
% vector) and .point (3x1, a point the plane passes through), i.e. the
% plane {r : dot(r - point, normal) = 0} -- plot_cross_section.m's
% general (non-axis-name) calling form.
%
%   latPlane.normal = [0 0 1]'                        (horizontal)
%   latPlane.point  = [0 0 refPoint(3)]'              (at ref's height)
%   lonPlane.normal = [-sin(phi) cos(phi) 0]'         (a vertical meridian)
%   lonPlane.point  = [0 0 0]'                        (always through the centre)
% where phi = atan2(refPoint(2), refPoint(1)) is the reference's own
% azimuthal angle about the polar (z) axis.

if isscalar(ref)
    ic = ref;
    if ic < 1 || ic > S.n_cell || ~S.cell_alive(ic)
        error('ring_planes_for_location:badcell', ...
              'Cell %d is not a valid alive cell in this snapshot.', ic);
    end
    n = S.cell_n(ic);
    vids = S.cell_vlist(ic, 1:n);
    refPoint = mean(S.r_api(vids, :), 1)';
else
    refPoint = ref(:);
    if numel(refPoint) ~= 3
        error('ring_planes_for_location:badref', ...
              'ref must be a single cell index or a 3-element [x y z] point.');
    end
end

phi = atan2(refPoint(2), refPoint(1));

latPlane.normal = [0; 0; 1];
latPlane.point  = [0; 0; refPoint(3)];

lonPlane.normal = [-sin(phi); cos(phi); 0];
lonPlane.point  = [0; 0; 0];
end
