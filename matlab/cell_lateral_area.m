function lat_area = cell_lateral_area(S)
% CELL_LATERAL_AREA  Total lateral (cell-cell wall) surface area for
% every ALIVE cell in S, computed directly from vertex geometry -- the
% snapshot format (read_vertex_snapshot.m) stores apical/basal area but
% not this, so it has to be built here from r_api/r_bas/cell_vlist.
%
% Same triangulation convention as Force.f90's own lateral-face term
% (and plot_cross_section.m's rendering of it): each lateral quad face
% (Ak,Bk,Bkp,Akp) split into two triangles, (Ak,Bk,Bkp) and
% (Ak,Bkp,Akp), summed over all n edges of the cell's own ring.
%
%   lat_area = cell_lateral_area(S)
%       one value per alive cell, in the same order as
%       find(S.cell_alive) -- i.e. matching tissue_color_values.m's own
%       per-cell ordering, so the two can be combined directly.

idx = find(S.cell_alive);
lat_area = zeros(numel(idx), 1);
for ii = 1:numel(idx)
    ic = idx(ii);
    n = S.cell_n(ic);
    vids = S.cell_vlist(ic, 1:n);
    Aapi = S.r_api(vids, :);
    Abas = S.r_bas(vids, :);
    a = 0.0;
    for k = 1:n
        kp = mod(k, n) + 1;
        a = a + 0.5 * norm(cross(Abas(k, :)  - Aapi(k, :), Abas(kp, :) - Aapi(k, :)));
        a = a + 0.5 * norm(cross(Abas(kp, :) - Aapi(k, :), Aapi(kp, :) - Aapi(k, :)));
    end
    lat_area(ii) = a;
end
end
