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

The example writes `examples/hall_channel.vtu`, readable by ParaView. Run the
test suite with:

```sh
julia --project=. test/runtests.jl
```

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
options = SolverOptions(max_area=2.5e-5, azimuthal_divisions=64)
result = solve(geometry, boundaries, gas; options)
write_vtk("hall_channel.vtu", result)
```

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
boundary segment is revolved into triangular surface patches. Direct solid
angles use the Oosterom–Strackee formula; a bounding-volume hierarchy removes
occluded rays. Axisymmetry reduces boundary exchange integration to one
receiver wedge and a full ring of emitters.

The segment exchange matrix is symmetrized for reciprocity and balanced to
satisfy enclosure closure. Prescribed source fluxes drive a linear diffuse
radiosity solve. Visible source distributions then give cell-centered number
density and the exact vector solid-angle moment gives mean velocity. The
`radiosity_residual`, `particle_balance_residual`, and
`exchange_closure_error` fields on `FlowResult` should be inspected for every
new geometry.

Increasing `azimuthal_divisions` resolves the surfaces of revolution more
accurately. Decreasing `max_area` refines the R–Z output mesh and boundary
segmentation. A zero `max_area` selects one two-hundredth of the polygon area.

## VTK fields

The `.vtu` file contains triangle cell data:

- `number_density` in m⁻³;
- `velocity` in m/s, with components `(u_z,u_r,u_θ)`;
- `direct_view_factor_<label>`, the unobstructed solid-angle fraction
  \(\Omega/(4\pi)\) for each physical boundary label;
- `density_from_<label>` in m⁻³, the final density emitted by that label after
  diffuse-wall coupling.

Field-data scalars record the three solver diagnostics. Labels are lowercased
and sanitized for use as VTK array names.

## Scope and physical assumptions

This release assumes neutral-neutral mean free paths are much longer than the
device dimensions. It is steady and collisionless and does not include Katz’s
time-dependent first-order transport equation or an ionization sink. It also
does not model drifting reservoirs, specular or partially accommodating walls,
multiple species, or neutral-neutral collisions. A diffuse-only closed domain
has no unique absolute steady density and is rejected; include a finite-area
inlet or reservoir/opening.
