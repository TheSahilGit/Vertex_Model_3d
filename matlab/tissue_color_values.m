function [cval, clabel] = tissue_color_values(S, colorby)
% TISSUE_COLOR_VALUES  Per-cell scalar values (and an axis label) for
% the requested colouring mode, shared by plot_tissue_3d.m/
% plot_cross_section.m and make_tissue_video.m/
% make_cross_section_video.m so a single snapshot and a whole video use
% exactly the same definition.
%
%   [cval, clabel] = tissue_color_values(S, 'nsides')        % number of sides
%   [cval, clabel] = tissue_color_values(S, 'volume')        % V/V0 - 1        (strain)
%   [cval, clabel] = tissue_color_values(S, 'volume_abs')    % V               (absolute)
%   [cval, clabel] = tissue_color_values(S, 'area')          % Aapi/A0 - 1     (apical strain)
%   [cval, clabel] = tissue_color_values(S, 'apical_area')   % Aapi            (absolute)
%   [cval, clabel] = tissue_color_values(S, 'basal_area')    % Abas/A0bas - 1  (basal strain)
%   [cval, clabel] = tissue_color_values(S, 'basal_area_abs')% Abas            (absolute)
%   [cval, clabel] = tissue_color_values(S, 'lateral_area')  % lateral wall area (absolute)
%   [cval, clabel] = tissue_color_values(S, 'total_area')    % Aapi+Abas+Alat  (absolute)
%   [cval, clabel] = tissue_color_values(S, 'shapefactor')   % total_area / V^(2/3)
%
% 'lateral_area', 'total_area' and 'shapefactor' need the lateral wall
% area, which (unlike apical/basal area) is not stored in the snapshot
% itself, so it is computed directly from vertex geometry
% (cell_lateral_area.m) -- the modes noticeably more expensive than a
% plain lookup of an already-stored field; see make_tissue_video.m/
% make_cross_section_video.m for how a whole video avoids paying for
% this twice per frame. Lateral area has no natural "strain" version
% (Lambda_line has no rest area the way K_A/K_A_bas have A0/A0_bas), so
% it and the two quantities built from it are absolute only.
%
% 'shapefactor' is the 3-D analogue of the 2-D vertex model's
% dimensionless shape index p = P/sqrt(A) (Bi & Manning and others):
% here q = S_total/V^(2/3), minimised by a sphere: rounder/more compact
% cells have a lower q, more elongated/irregular ones a higher q.
%
% cval is one value per ALIVE cell, in the same order as
% find(S.cell_alive) (i.e. matching cell_faces_matrix.m's row order).

idx = find(S.cell_alive);
switch lower(colorby)
    case 'nsides'
        cval = S.cell_n(idx);
        clabel = 'number of sides';
    case 'volume'
        cval = S.cell_Vlast(idx) ./ max(S.cell_V0(idx), eps) - 1;
        clabel = 'volume strain  (V/V_0 - 1)';
    case 'volume_abs'
        cval = S.cell_Vlast(idx);
        clabel = 'cell volume';
    case 'area'
        cval = S.cell_Alast(idx) ./ max(S.cell_A0(idx), eps) - 1;
        clabel = 'apical area strain  (A/A_0 - 1)';
    case 'apical_area'
        cval = S.cell_Alast(idx);
        clabel = 'apical area';
    case 'basal_area'
        cval = S.cell_Abaslast(idx) ./ max(S.cell_A0bas(idx), eps) - 1;
        clabel = 'basal area strain  (A^{bas}/A_0^{bas} - 1)';
    case 'basal_area_abs'
        cval = S.cell_Abaslast(idx);
        clabel = 'basal area';
    case 'lateral_area'
        cval = cell_lateral_area(S);
        clabel = 'lateral (cell-cell) wall area';
    case 'total_area'
        cval = S.cell_Alast(idx) + S.cell_Abaslast(idx) + cell_lateral_area(S);
        clabel = 'total surface area  (apical + basal + lateral)';
    case 'shapefactor'
        Stot = S.cell_Alast(idx) + S.cell_Abaslast(idx) + cell_lateral_area(S);
        cval = Stot ./ max(S.cell_Vlast(idx), eps).^(2/3);
        clabel = 'shape factor  (S_{total} / V^{2/3})';
    otherwise
        error('tissue_color_values:colorby', 'unknown colorby option "%s"', colorby);
end
end
