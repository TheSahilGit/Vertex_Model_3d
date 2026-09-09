function [cval, clabel] = tissue_color_values(S, colorby)
% TISSUE_COLOR_VALUES  Per-cell scalar values (and an axis label) for
% the requested colouring mode, shared by plot_tissue_3d.m and
% make_tissue_video.m so a single snapshot and a whole video use
% exactly the same definition.
%
%   [cval, clabel] = tissue_color_values(S, 'nsides')
%   [cval, clabel] = tissue_color_values(S, 'volume')
%   [cval, clabel] = tissue_color_values(S, 'area')
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
    case 'area'
        cval = S.cell_Alast(idx) ./ max(S.cell_A0(idx), eps) - 1;
        clabel = 'apical area strain  (A/A_0 - 1)';
    otherwise
        error('tissue_color_values:colorby', 'unknown colorby option "%s"', colorby);
end
end
