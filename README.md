# 3D Vertex Model on a Hollow Spherical Shell

A Fortran 90 implementation of a 3D vertex model for an epithelial-like
tissue tiling the surface of a sphere, with Matlab visualization.

## The model

**Geometry.** Cells are polyhedra sandwiched between two concentric
spherical shells: an outer **apical** surface (radius `R_apical`) and
an inner **basal** surface (radius `R_basal`). The shell is **hollow**
— there is nothing modelled at `r < R_basal`, only the layer of cells
between the two radii. Every "vertex column" `j` carries one apical
point `r_api(:,j)` and one basal point `r_bas(:,j)` directly below it;
a cell is the polyhedron whose apical face is the polygon through its
apical points, whose basal face is the corresponding basal polygon,
and whose lateral faces are the quadrilaterals joining consecutive
apical/basal edges. This is the standard columnar-epithelium
simplification of the 3D vertex model (Honda/Okuda-type).

**Initial mesh.** The initial tissue is a geodesic (subdivided
icosahedral) triangulation of the sphere, whose *dual* tessellation
gives the cells: one cell per triangulation vertex, one vertex column
per triangulation face. This guarantees exactly 12 pentagonal cells
and otherwise hexagonal cells, i.e. a proper "buckyball-like" packing,
and every interior vertex column is shared by exactly 3 cells
(trivalent), as required by the vertex-model formalism.

**Energy functional** (standard terms):

```
E = Σ_cells (K_V/2)(V     - V0)²          volume elasticity
  + Σ_cells (K_A/2)(Aapi  - A0)²          apical-area elasticity
  + Σ_cells (K_P/2) Papi²                apical perimeter contractility
  + Σ_cells (K_A_bas/2)(Abas - A0_bas)²  basal-area elasticity
  + Σ_cells (K_P_bas/2) Pbas²            basal perimeter contractility
  + Σ_edges Λ · S_lateral(edge)          lateral (cell-cell) interfacial tension
```

The basal-area/perimeter terms are the exact mirror of the apical
ones, built from the basal shell instead of the apical one. They are
genuinely independent moduli (`K_A_bas`, `K_P_bas`, and a separate
target-area scale `A0_bas_scale`) — set them equal to their apical
counterparts for a symmetric apical/basal tissue, or differently for
an asymmetric one (e.g. a stiffer or more contractile basal side). If
omitted from `para.in`, each defaults to its apical counterpart
(`K_A_bas=K_A`, `K_P_bas=K_P`, `A0_bas_scale=A0_scale`), so older
`para.in` files still work unchanged, with the basal face behaving
exactly like the apical one.

Volumes and areas are computed by triangulating every face (fan from
the face centroid, or two triangles per lateral quad) and using exact
divergence-theorem/triangle-area identities; every term is
differentiated **analytically** (see `src/Force.f90`), not by finite
differences. `src/sanity_checks.f90` verifies the analytic gradients
against finite differences and verifies the whole volume/area pipeline
(now including the basal area) against an exact unit cube.

**Dynamics.** Overdamped Langevin: `ζ dr/dt = -dE/dr + ξ(t)`, Euler–
Maruyama integrated (`src/Langevin_update.f90`).

**Topology (junctional rearrangements)**, all in terms of the ring
data structure `cells(:)%vlist` (`src/mod_topology.f90`):
- **T1** (`src/T1_transition.f90`): edges shorter than
  `L_T1_threshold` are flipped between the four surrounding cells.
- **T2** (`src/T2_transition.f90`): triangular cells whose area drops
  below `A_T2_threshold` collapse to a point and are removed; this can
  cascade (a collapse can leave a neighbour with < 3 sides), which is
  handled by a bounded queue of forced follow-up collapses.
- **Cell division** (`src/Cell_Division.f90`): a cell whose volume
  exceeds `V_division_threshold · V0` splits into two daughters along
  a roughly-opposite pair of edges.
- **T4** (`src/T4_transition.f90`, optional, off by default): crowding-
  induced live-cell extrusion — a cell squeezed below
  `V_extrusion_threshold · V0` is removed from the tissue, but (unlike
  T2) not by collapsing to a point: it is whittled down one ordinary T1
  flip of its own shortest edge per check (so its neighbours merge
  faces around it, gaining new adjacencies one at a time), only reusing
  T2's point-collapse once it has been reduced to a triangle. Purely
  reuses the existing T1/T2 primitives — no new ring-surgery code.
  Extrusion direction (apical- vs basal-ward) is classified per event
  by comparing apical vs basal compression (`A_last/A0` vs
  `A_bas_last/A0_bas`) at the moment of removal.

All four are implemented so that every mutation is validated (edge
count, capacity) **before** any global state is touched, so a failed
attempt never leaves a half-applied/corrupted mesh. Correctness is
checked continuously at runtime: an Euler-characteristic check
(`V - E + F = 2`) and a full ring-integrity check (no dangling/dead/
duplicate vertex references) run after every batch of topology events
and division events, aborting the run if either ever fails.

## Layout

```
para.in               simulation parameters (Fortran NAMELIST, heavily commented)
src/                  Fortran 90 source (see Makefile for the compile order)
  mod_kinds.f90         precision kinds
  mod_parameters.f90    parameter file I/O
  mod_data.f90          core data structures (cells, vertex arrays)
  mod_geometry.f90      triangle-area / tetrahedron-volume + analytic gradients
  mod_topology.f90      ring surgery + Euler-characteristic check
  mesh_init.f90         icosahedral geodesic mesh + dual tessellation
  Force.f90             energy functional + forces
  Langevin_update.f90   overdamped Langevin integrator
  T1_transition.f90     T1 neighbour exchange
  T2_transition.f90     T2 cell removal (with cascade handling)
  T4_transition.f90     T4 crowding-induced extrusion (optional, off by default)
  Cell_Division.f90     cell division
  io_module.f90         binary snapshot writer (data/)
  sanity_checks.f90     start-up self-tests + runtime invariant checks
  main.f90              driver
build/                object files + .mod files (created by `make`)
data/                 binary snapshots + metadata (created at runtime)
matlab/               Matlab analysis & plotting (see below)
```

## Building and running

```sh
make          # builds ./vertex3d
./vertex3d    # reads para.in, writes data/*.dat
```

`make clean` removes `build/` and the executable; `make distclean` also
removes `data/`.

## Output format

Each `data/snap_<step>.dat` is a flat, **stream-access** unformatted
binary file (no Fortran record markers), documented byte-for-byte in
the header comment of `src/io_module.f90`. `data/mesh_meta.txt` is a
small plain-text key/value metadata file, and `data/dump_list.txt`
lists every snapshot written, in order.

Alongside every `snap_<step>.dat`, a small `data/diag_<step>.dat` (same
`<step>` and cadence, so the two always pair up) holds scalar
whole-tissue diagnostics for that step: total energy, lumen volume
(enclosed by the basal/inner surface), outer surface area, the largest
`|force|` over every vertex, the alive cell count, the cumulative
T1/T2 event counts since t=0, and (if T4 is enabled) the cumulative T4
extrusion count and its apical-ward/basal-ward/ambiguous breakdown.
`data/diag_list.txt` lists them the same way `dump_list.txt` does for
snapshots.

`data/` itself is regenerated by every run and gitignored, except
`data/example/` — a small bundled example kept under version control
so the Matlab scripts below have something to run against without
building and running the Fortran code first: three `snap_`/`diag_`
pairs (first/middle/last step) for the videos, plus the *entire*
`diag_` series from that same short run (each file is under 100 bytes,
so keeping all of them costs nothing) for a properly-populated example
time-series plot.

## Matlab analysis (`matlab/`)

There are two independent driver scripts:

- **`main_analysis.m`** reads the scalar `diag_<step>.dat` time series
  and plots total energy, lumen volume, outer surface area, max force,
  cell count, cumulative T1/T2 counts, and cumulative T4 (with its
  apical/basal/ambiguous breakdown), all against time, in one
  multi-panel figure.
- **`main_movie.m`** reads the vertex `snap_<step>.dat` series and
  renders up to THREE independent videos, each its own figure/window:
  the whole tissue, a latitude ring, and a longitude ring through a
  chosen cell/location. The three are fully separate (no shared
  figure), so each can be freely rotated/inspected on its own without
  affecting the others.

```matlab
main_analysis.m            % driver: diagnostics-vs-time plot (see below)
main_movie.m               % driver: up to three independent videos (see below)
read_vertex_snapshot.m     % read one binary vertex snapshot -> struct
read_diagnostics.m         % read one binary diagnostics snapshot -> struct
read_mesh_meta.m           % read data/mesh_meta.txt -> struct
list_snapshots.m           % list available vertex snapshots, with their
                            % saved iteration numbers, in order
list_diagnostics.m         % (same, for the diagnostics series)
select_snapshots.m         % pick out an array (or scalar) of saved
                            % iteration numbers from a list_snapshots.m result
cell_faces_matrix.m        % snapshot -> patch()-ready Faces matrix
cell_lateral_area.m        % snapshot -> per-cell lateral (cell-cell wall) area,
                            % computed from raw geometry (not stored in the
                            % snapshot the way apical/basal area are)
tissue_color_values.m      % snapshot + colouring mode -> per-cell scalars --
                            % 'nsides', 'volume'/'volume_abs', 'area'/
                            % 'apical_area', 'basal_area'/'basal_area_abs',
                            % 'lateral_area', 'total_area', 'shapefactor'
                            % (the 3-D shape index S_total/V^(2/3))
plot_tissue_3d.m           % standard whole-tissue 3D plot (colorbar, title)
plot_cross_section.m       % ring cross-section: hollow interior + individual
                            % cell shapes, per-cell coloured; cut by either an
                            % axis-aligned plane (legacy 'x'|'y'|'z' + value)
                            % or an arbitrary plane (normal + point struct)
ring_planes_for_location.m % a cell index or [x y z] point -> the latitude
                            % and longitude cutting planes through it
make_tissue_video.m        % render a whole-tissue video across many snapshots
make_cross_section_video.m % (same, for a single cross-section view --
                            % called once each for the latitude/longitude rings)
```

```matlab
S = read_vertex_snapshot('data/example/snap_00005000.dat');
plot_tissue_3d(S, 'nsides');                     % or 'volume', 'shapefactor', ...
plot_tissue_3d(S, 'shapefactor', 'ColorMap', 'turbo');  % any built-in colormap name
plot_cross_section(S, 'y', 0.0);   % legacy form: cut through the sphere centre

[latPlane, lonPlane] = ring_planes_for_location(S, 7);   % cell 7's own location
plot_cross_section(S, latPlane, 'ColorBy', 'volume', 'Title', 'Latitude ring');
plot_cross_section(S, lonPlane, 'ColorBy', 'volume', 'Title', 'Longitude ring');

D = read_diagnostics('data/example/diag_00005000.dat');
D.energy, D.lumen_volume, D.outer_area, D.max_force, D.n_cells, ...
D.cumulative_T1, D.cumulative_T2, D.cumulative_T4, ...
D.cumulative_T4_apical, D.cumulative_T4_basal, D.cumulative_T4_ambiguous
```

Both driver scripts locate their own folder (work whether you `cd`
into `matlab/` and run them directly, or run them from the project
root via `run('matlab/main_analysis.m')` / `run('matlab/main_movie.m')`),
read from `data/`, and automatically fall back to the bundled
`data/example/` if `data/` is empty — so on a fresh clone, before ever
building or running the Fortran code, either one alone already
produces something:

```matlab
% from the project root
run('matlab/main_analysis.m')   % diagnostics-vs-time figure
run('matlab/main_movie.m')      % up to three independent videos
```

A block of flags at the top of `main_movie.m` controls what it makes.
`flag_tissue`/`flag_latitude`/`flag_longitude` are fully independent —
run any subset, in any combination:

```matlab
flag_tissue     = true;    % make the whole-tissue (outside) view
flag_latitude   = true;    % make the latitude-ring cross-section view
flag_longitude  = true;    % make the longitude-ring cross-section view

colorby = 'nsides';    % colouring for all three -- see tissue_color_values.m
                        % for the full list ('nsides', 'volume'/'volume_abs',
                        % 'area'/'apical_area', 'basal_area'/'basal_area_abs',
                        % 'lateral_area', 'total_area', 'shapefactor')
colormapName = 'parula'; % any built-in MATLAB colormap name, e.g. 'turbo',
                        % 'jet', 'hot', 'cool', 'copper', 'bone'
ref_location = 1;       % which cell/location the latitude and longitude
                        % rings are cut through -- EITHER a cell index
                        % (e.g. 1) OR an explicit [x y z] point (e.g.
                        % [5 0 9]). Resolved once, from the FIRST
                        % selected frame, into two fixed cutting planes
                        % reused for every later frame.
video_its = [];         % which saved snapshots -- see below
video_fps = 8;
```

"Latitude" and "longitude" are defined relative to the z-axis as the
polar axis (this mesh has no intrinsically marked pole; z is simply
the axis every other plot here already treats as "up") --
see `ring_planes_for_location.m`: the latitude ring is the horizontal
plane at the reference location's own height, the longitude ring is
the vertical meridian plane through the poles and that location's own
azimuthal angle. A reference point very close to a pole can end up
with no cells straddling one of these planes at all (an empty plot,
with a printed warning); and at coarse mesh resolutions, a latitude
ring can legitimately come out with visible gaps between clusters of
cells rather than one continuous loop -- a real property of exactly
where that specific cut falls on a coarse icosahedral mesh, not a bug
(the longitude ring through the same cell, or a finer mesh, is usually
a clean unbroken loop).

`video_its` is an array given directly in terms of the **saved
iteration numbers** (the `it_dumps` cadence from `para.in`,
not raw simulation steps or a plain file-list index) — e.g.
`video_its = 1000:100:5000`. A scalar, e.g. `video_its = 5000`, works
exactly like a 1-element array: the same code path still "makes the
video", which in that case is just that one plot. Leave it `[]` to use
every available snapshot. A requested iteration that wasn't actually
saved is matched to the nearest one that was, with a printed note.

Videos are written to `videos/` (gitignored) as `.avi` (`Motion JPEG
AVI`, chosen because it — unlike `MPEG-4` — is supported by MATLAB on
every platform, including Linux). `make_tissue_video.m`/
`make_cross_section_video.m` reuse one figure for their whole run
(left open on the last frame so you can look at it), pin it to an
exact pixel *size* every frame (never its on-screen position, which is
left free to drag around) and force-resize any captured frame that
still comes back a different size before handing it to `VideoWriter`,
since `writeVideo` errors out on the first size mismatch (e.g. a
colorbar tick label gaining a digit can shift the rendered axes by a
pixel) — otherwise a good chunk of frames into a long video. With more
than one flag on, the videos are made one after another, not at the
same time: this MATLAB install is licensed for Parallel Computing
Toolbox but does not actually have it installed (checked directly --
`parpool`/`gcp` are undefined here), so there is no way to run them
concurrently within one MATLAB session.

## Sanity checks built in

- Analytic vs. finite-difference gradients for both the triangle-area
  and tetrahedron-volume primitives (`sanity_checks.f90`).
- A unit-cube volume/area test run through the *production*
  `compute_forces` routine, to pin down the face-winding convention.
- Geodesic-mesh vertex/face counts and pentagon/hexagon degree
  distribution checked exactly against the closed-form icosahedral
  formulas.
- Euler characteristic (`V - E + F = 2`) and full ring integrity
  (no dead/duplicate/dangling vertex references) checked after every
  batch of T1/T2/division events and at the end of the run.
- Total apical area and total shell volume compared against the ideal
  continuum sphere/shell values at start-up.

## Known limitations

- Polygon valence is capped at `MAX_SIDES = 16` (`src/mod_data.f90`);
  an event that would exceed this is skipped with a warning.
- Cell/vertex array capacity is pre-allocated at
  `capacity_growth_factor × (initial count)`; exceeding it likewise
  skips further divisions with a warning (raise the parameter if you
  see this often).
- T2 cascade handling is bounded (`QCAP` in `T2_transition.f90`) for
  safety on pathological meshes; ordinary runs never approach this.
