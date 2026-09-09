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
E = Σ_cells (K_V/2)(V   - V0)²      volume elasticity
  + Σ_cells (K_A/2)(Aapi - A0)²     apical-area elasticity
  + Σ_cells (K_P/2) Papi²           apical perimeter contractility
  + Σ_edges Λ · S_lateral(edge)     lateral (cell-cell) interfacial tension
```

Volumes and areas are computed by triangulating every face (fan from
the face centroid, or two triangles per lateral quad) and using exact
divergence-theorem/triangle-area identities; every term is
differentiated **analytically** (see `src/Force.f90`), not by finite
differences. `src/sanity_checks.f90` verifies the analytic gradients
against finite differences and verifies the whole volume/area pipeline
against an exact unit cube.

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

All three are implemented so that every mutation is validated (edge
count, capacity) **before** any global state is touched, so a failed
attempt never leaves a half-applied/corrupted mesh. Correctness is
checked continuously at runtime: an Euler-characteristic check
(`V - E + F = 2`) and a full ring-integrity check (no dangling/dead/
duplicate vertex references) run after every batch of topology events
and division events, aborting the run if either ever fails.

## Layout

```
para_Simulation.dat   simulation parameters (Fortran NAMELIST, heavily commented)
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
./vertex3d    # reads para_Simulation.dat, writes data/*.dat
```

`make clean` removes `build/` and the executable; `make distclean` also
removes `data/`.

## Output format

Each `data/snap_<step>.dat` is a flat, **stream-access** unformatted
binary file (no Fortran record markers), documented byte-for-byte in
the header comment of `src/io_module.f90`. `data/mesh_meta.txt` is a
small plain-text key/value metadata file, and `data/dump_list.txt`
lists every snapshot written, in order.

`data/` itself is regenerated by every run and gitignored, except
`data/example/` — three snapshots (first/middle/last step) from a
default-parameter run, kept under version control so the Matlab
scripts below have something to run against without building and
running the Fortran code first.

## Matlab analysis (`matlab/`)

```matlab
main_analysis.m          % example driver: reads the latest snapshot,
                          % produces the whole-tissue plot and the
                          % cutaway cross-section
read_vertex_snapshot.m   % read one binary snapshot -> struct
read_mesh_meta.m         % read data/mesh_meta.txt -> struct
list_snapshots.m         % list available snapshots in order
cell_faces_matrix.m      % snapshot -> patch()-ready Faces matrix
plot_tissue_3d.m         % standard whole-tissue 3D plot
plot_cross_section.m     % cutaway view: hollow interior + ring +
                          % individual cell shapes
```

```matlab
S = read_vertex_snapshot('data/example/snap_00005000.dat');
plot_tissue_3d(S, 'nsides');       % or 'volume' / 'area'
plot_cross_section(S, 'y', 0.0);   % cut through the sphere centre
```

`main_analysis.m` locates its own folder (works whether you `cd` into
`matlab/` and run it directly, or run it from the project root via
`run('matlab/main_analysis.m')`), reads the latest snapshot from
`data/`, and automatically falls back to the bundled
`data/example/` if `data/` is empty — so on a fresh clone, before
ever building or running the Fortran code, this alone already
produces all four plots:

```matlab
% from the project root
run('matlab/main_analysis.m')
```

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
