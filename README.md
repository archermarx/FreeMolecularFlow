# FreeMolecularFlow.jl

`FreeMolecularFlow.jl` is a deterministic, axisymmetric neutral-gas solver for
the free-molecular regime of plasma thrusters. It implements a steady,
collisionless specialization of the view-factor algorithm introduced by Katz
and Mikellides:

- I. Katz and I. G. Mikellides, [“A New Algorithm for the Neutral Gas in the
  Free-Molecule Regimes of Hall and Ion Thrusters,” IEPC-2009-095](https://electricrocket.org/IEPC/IEPC-2009-095.pdf).
- I. Katz and I. G. Mikellides, [“Neutral gas free molecular flow algorithm
  including ionization and walls for use in plasma simulations,” JCP 230,
  1454–1464 (2011)](https://doi.org/10.1016/j.jcp.2010.11.013).

The code generates no Monte Carlo particle noise. It supports a single neutral
species, prescribed mass-flow inlets, stationary back-pressure reservoirs,
vacuum outlets, and fully accommodating diffuse walls.

## Install and run

Use Julia 1.10 or newer:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. bin/free_molecular_flow.jl examples/hall_channel.toml
```

Longer cases can print periodic progress and residuals by passing an interval
in seconds:

```sh
julia --project=. bin/free_molecular_flow.jl --status-interval 5 examples/spt100.toml
```

The compact status table reports the active phase, phase-local
iteration/progress, total elapsed time, and every residual available at that
point. The default interval is zero, which disables status output. The library equivalents are
`solve(...; status_interval=5)` and `run_config(...; status_interval=5)`.

The example writes `examples/hall_channel.vtu`, readable by ParaView. Run the
test suite with:

```sh
julia --project=. test/runtests.jl
```

For a runnable input file with comments describing the main geometry,
boundary, solver, cache, output, and line-extraction options, see
[`examples/annotated_options.toml`](examples/annotated_options.toml).

## Julia API

Geometry is one simple closed polygon in the `(z,r)` plane. Edge `i` joins
point `i` to point `i+1` cyclically, and has the corresponding label. All input
and output use SI units.

```julia
using FreeMolecularFlow

geometry = AxisymmetricGeometry(
    [(0.0,0.0), (0.05,0.0), (0.05,0.03), (0.0,0.03)],
    ["axis", "vacuum", "wall", "anode"])

boundaries = Dict(
    "axis"   => Axis(),
    "vacuum" => BackPressure(0.0, 300.0),
    "wall"   => DiffuseWall(500.0),
    "anode"  => Inflow(5e-6, 700.0))

gas = Gas(131.293; unit=:amu)
options = SolverOptions(max_area=2.5e-5, max_boundary_length=5e-3,
                        azimuthal_divisions=64)
result = solve(geometry, boundaries, gas; options)
write_vtk("hall_channel.vtu", result)
```

For parameter sweeps, prepare the geometry-dependent transport operator once:

```julia
prepared = prepare(geometry, boundaries; options)
result_1 = solve(prepared, gas)

updated = copy(boundaries)
updated["anode"] = Inflow(7e-6, 700.0)
result_2 = solve(prepared, updated, gas)
```

Prepared data can also be reused across Julia processes:

```julia
prepared = prepare_cached("spt100.fmf-cache",geometry,boundaries;options)
result = solve(prepared,boundaries,gas)
```

or directly from TOML:

```toml
[cache]
path = "spt100.fmf-cache"
```

Relative cache paths are case-relative. A missing cache is built atomically;
one whose geometry or solver options changed is rebuilt. Cache files are tied
to the Julia and FreeMolecularFlow versions that wrote them, because Julia's
native serialization format is not intended as a portable interchange format.

The mesh, revolved surface, visibility hierarchy, boundary-exchange matrix, and
cell solid-angle moments are reused. Replacement boundaries must use the same
geometry labels, but their values and types may change.

Polygon orientation is normalized internally. The polygon must not intersect
itself, all radii must be nonnegative, and any edge wholly on `r=0` must be an
`Axis()` edge. The current release deliberately supports one loop only; holes
and disconnected flow regions are rejected by construction.

### Boundary conditions

- `Inflow(mass_flow_rate, temperature)` distributes the total kg/s uniformly
  over every revolved segment with the same label and emits a cosine-law
  half-Maxwellian.
- `BackPressure(pressure, temperature)` admits the inward half of a stationary
  Maxwellian reservoir and permits incident particles to escape. Use zero Pa
  for a vacuum opening.
- `DiffuseWall(temperature)` completely accommodates incident particles and
  diffusely re-emits the same number flux at the wall temperature.
- `Axis()` represents symmetry, not a physical surface, and has zero area.

Specular reflection is not implemented. A TOML boundary with a specular type
fails with an explicit diagnostic rather than silently changing its physics.

## Numerical method

The polygon is refined to a constrained Delaunay triangle mesh. Every physical
boundary segment is revolved into triangular solid-angle integration patches.
Direct solid angles use the Oosterom–Strackee formula. Visibility rays are
intersected analytically with the original edges revolved into finite cylinders,
annular disks, and cones; a small hierarchy over those exact surfaces rejects
distant geometry. Axisymmetry reduces boundary exchange integration to one
receiver wedge and a full ring of emitters. Visibility cost therefore scales
with the number of input polygon edges rather than the azimuthal patch count.

The segment exchange matrix is symmetrized for reciprocity and balanced to
satisfy enclosure closure. Prescribed source fluxes drive a linear diffuse
radiosity solve. Visible source distributions then give cell-centered number
density, mean velocity, and centered directional temperatures. Exact scalar,
vector, and diagonal second solid-angle moments are evaluated over each visible
spherical triangle, so temperature does not assume a local Maxwellian VDF. The
`radiosity_residual`, `particle_balance_residual`, and
`exchange_closure_error` fields on `FlowResult` should be inspected for every
new geometry.

Increasing `azimuthal_divisions` resolves the surfaces of revolution more
accurately. Decreasing `max_area` refines the R–Z output mesh. A zero `max_area`
selects one two-hundredth of the polygon area.
`max_boundary_length` independently limits the R–Z length of transport boundary
segments. Set it explicitly when refining the volume mesh so the exchange
matrix and revolved source count remain fixed; zero selects an automatic value
of twice the characteristic cell length.
Even values of `azimuthal_divisions` are preferred: their mirrored surface
quadrature lets field reconstruction evaluate half of each ring. Odd values
remain supported and use the full ring.

Azimuthal resolution can also be selected automatically. Set
`azimuthal_tolerance` to a positive relative tolerance and use
`azimuthal_divisions` as the starting resolution. The solver doubles that
resolution until the exchange matrix and all scalar, first, and second
directional moments satisfy the tolerance, or throws when
`max_azimuthal_divisions` is reached. Both limits must be even in adaptive
mode. For example:

```toml
[solver]
azimuthal_divisions = 16
azimuthal_tolerance = 5.0e-2
max_azimuthal_divisions = 128
```

The actual resolution and estimated change are available as
`result.azimuthal_divisions` and `result.azimuthal_convergence_error`, and are
also written as VTK field data. The default tolerance is zero, retaining fixed
resolution and avoiding the cost of convergence comparisons.

Field reconstruction uses every Julia worker thread made available at process
startup. For example, run the command-line solver on the automatically selected
thread count with:

```sh
julia --threads=auto --project=. bin/free_molecular_flow.jl case.toml
```

## VTK fields

The `.vtu` file contains triangle cell data:

- `number_density` in m⁻³;
- `velocity` in m/s, with components `(u_z,u_r,u_θ)`;
- `temperature` in K, with centered second-moment components `(T_z,T_r,T_θ)`;
- `direct_view_factor_<label>`, the unobstructed solid-angle fraction
  \(\Omega/(4\pi)\) for each physical boundary label;
- `density_from_<label>` in m⁻³, the final density emitted by that label after
  diffuse-wall coupling.

Field-data scalars record the three conservation/residual diagnostics plus the
actual azimuthal resolution and its convergence estimate. Labels are
lowercased and sanitized for use as VTK array names.

## Line extraction

Add one or more array-of-table entries to a TOML case. Two control points make
a straight line; additional points make a polyline:

```toml
[[extraction_lines]]
name = "plume_centerline"
points = [
  { z = 0.025, r = 0.0 },
  { z = 0.100, r = 0.0 },
]
num_points = 101
method = "direct"
outside_domain = "keep"
fields = ["number_density", "velocity", "temperature", "view_factors", "density_contributions"]
filename = "spt100_centerline.csv"
```

Sampling is uniform in cumulative path distance and includes both endpoints.
Specify exactly one of `num_points` or `spacing` (metres). With spacing, a final
short interval is added as needed so the last control point is always included.

`method = "direct"` (the default) evaluates the view-factor moments at the exact
sample coordinates. `method = "cell"` is faster and copies the containing
triangle's cell-centered result. `outside_domain` may be `"keep"`, `"drop"`, or
`"error"`; kept outside points have `cell_index=0` and empty solution fields.

The CSV always includes sample number, polyline segment, cumulative fraction
and distance (metres), `(z,r)` coordinates (metres), domain membership, and
containing cell index. The available solution-field groups are `number_density`,
`velocity`, `temperature`, `view_factors`, and `density_contributions`.
Omitting `fields` writes all groups. Omitting
`filename` derives `<name>.csv`; otherwise the path is resolved relative to the
TOML file. Multiple `[[extraction_lines]]` blocks produce multiple CSV files.
When CLI status output is enabled, direct extraction progress appears as an
`extract` phase in the same status table.

## Scope and physical assumptions

This release assumes neutral-neutral mean free paths are much longer than the
device dimensions. It is steady and collisionless and does not include Katz’s
time-dependent first-order transport equation or an ionization sink. It also
does not model drifting reservoirs, specular or partially accommodating walls,
multiple species, or neutral-neutral collisions. A diffuse-only closed domain
has no unique absolute steady density and is rejected; include a finite-area
inlet or reservoir/opening.
