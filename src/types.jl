const BOLTZMANN = 1.380649e-23
const ATOMIC_MASS = 1.66053906660e-27

abstract type BoundaryCondition end

"""A prescribed cosine-law inflow with total mass flow rate (kg/s)."""
struct Inflow <: BoundaryCondition
    mass_flow_rate::Float64
    temperature::Float64
    function Inflow(mass_flow_rate::Real, temperature::Real)
        mass_flow_rate >= 0 || throw(ArgumentError("inflow mass_flow_rate must be nonnegative"))
        temperature > 0 || throw(ArgumentError("inflow temperature must be positive"))
        new(Float64(mass_flow_rate), Float64(temperature))
    end
end

"""An open boundary connected to a stationary reservoir at pressure (Pa) and temperature (K)."""
struct BackPressure <: BoundaryCondition
    pressure::Float64
    temperature::Float64
    function BackPressure(pressure::Real, temperature::Real)
        pressure >= 0 || throw(ArgumentError("back-pressure pressure must be nonnegative"))
        temperature > 0 || throw(ArgumentError("back-pressure temperature must be positive"))
        new(Float64(pressure), Float64(temperature))
    end
end

"""A fully accommodating, cosine-law diffuse wall."""
struct DiffuseWall <: BoundaryCondition
    temperature::Float64
    function DiffuseWall(temperature::Real)
        temperature > 0 || throw(ArgumentError("wall temperature must be positive"))
        new(Float64(temperature))
    end
end

"""A zero-area symmetry-axis edge. Axis edges must lie at r = 0."""
struct Axis <: BoundaryCondition end

"""A single-species neutral gas. `molecular_mass` is in kg per molecule."""
struct Gas
    molecular_mass::Float64
    function Gas(molecular_mass::Real; unit::Symbol=:kg)
        mass = unit === :kg ? Float64(molecular_mass) :
               unit === :amu ? Float64(molecular_mass) * ATOMIC_MASS :
               throw(ArgumentError("Gas unit must be :kg or :amu"))
        mass > 0 || throw(ArgumentError("molecular mass must be positive"))
        new(mass)
    end
end

"""Mean molecular speed of a Maxwellian gas."""
mean_molecular_speed(gas::Gas, temperature::Real) =
    sqrt(8BOLTZMANN * Float64(temperature) / (pi * gas.molecular_mass))

"""
    AxisymmetricGeometry(points, edge_labels)

A simple closed polygon in the `(z,r)` meridional plane. `edge_labels[i]`
labels the edge from `points[i]` to `points[mod1(i+1,n)]`.
"""
struct AxisymmetricGeometry
    points::Vector{NTuple{2,Float64}}
    edge_labels::Vector{String}
    function AxisymmetricGeometry(points, edge_labels)
        pts = NTuple{2,Float64}[(Float64(p[1]), Float64(p[2])) for p in points]
        labels = String.(edge_labels)
        length(pts) >= 3 || throw(ArgumentError("geometry needs at least three points"))
        length(labels) == length(pts) ||
            throw(ArgumentError("edge_labels must have one entry per polygon edge"))
        all(p -> all(isfinite, p), pts) || throw(ArgumentError("geometry coordinates must be finite"))
        all(p -> p[2] >= 0, pts) || throw(ArgumentError("axisymmetric radius r must be nonnegative"))
        any(isempty, labels) && throw(ArgumentError("boundary labels may not be empty"))
        _validate_polygon!(pts)
        area = _signed_area(pts)
        if area < 0
            reverse!(pts)
            # Reversing vertices maps new edge i to old edge n-i.
            labels = [labels[mod1(length(labels) - i, length(labels))] for i in 1:length(labels)]
        end
        new(pts, labels)
    end
end

Base.@kwdef struct SolverOptions
    max_area::Float64 = 0.0
    min_angle::Float64 = 20.0
    azimuthal_divisions::Int = 64
    radiosity_tolerance::Float64 = 1e-10
    visibility_tolerance::Float64 = 1e-10
    max_mesh_points::Int = 200_000
end

struct BoundarySegment
    a::NTuple{2,Float64}
    b::NTuple{2,Float64}
    label::String
    condition::BoundaryCondition
    area::Float64
end

struct RZMesh
    points::Vector{NTuple{2,Float64}}
    cells::Vector{NTuple{3,Int}}
    centers::Vector{NTuple{2,Float64}}
    volumes::Vector{Float64}
    boundary_segments::Vector{BoundarySegment}
end

"""Results of a steady free-molecular calculation."""
struct FlowResult
    mesh::RZMesh
    labels::Vector{String}
    direct_view_factors::Dict{String,Vector{Float64}}
    density_contributions::Dict{String,Vector{Float64}}
    density::Vector{Float64}
    velocity::Matrix{Float64}       # 3 × ncells, ordered (z,r,azimuthal)
    boundary_flux::Vector{Float64}  # particles m^-2 s^-1
    radiosity_residual::Float64
    particle_balance_residual::Float64
    exchange_closure_error::Float64
end

number_density(result::FlowResult) = result.density

function _signed_area(p)
    0.5 * sum(p[i][1] * p[mod1(i+1, length(p))][2] -
              p[mod1(i+1, length(p))][1] * p[i][2] for i in eachindex(p))
end

function _orient(a, b, c)
    (b[1]-a[1]) * (c[2]-a[2]) - (b[2]-a[2]) * (c[1]-a[1])
end

function _on_segment(a, b, p; atol=1e-13)
    abs(_orient(a,b,p)) <= atol * max(1.0, hypot(b[1]-a[1], b[2]-a[2])) &&
        min(a[1],b[1])-atol <= p[1] <= max(a[1],b[1])+atol &&
        min(a[2],b[2])-atol <= p[2] <= max(a[2],b[2])+atol
end

function _segments_intersect(a,b,c,d)
    o1, o2, o3, o4 = _orient(a,b,c), _orient(a,b,d), _orient(c,d,a), _orient(c,d,b)
    ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
        ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0)) && return true
    return _on_segment(a,b,c) || _on_segment(a,b,d) || _on_segment(c,d,a) || _on_segment(c,d,b)
end

function _validate_polygon!(p)
    n = length(p)
    abs(_signed_area(p)) > eps(Float64) || throw(ArgumentError("polygon has zero area"))
    for i in 1:n
        a, b = p[i], p[mod1(i+1,n)]
        a != b || throw(ArgumentError("polygon has a duplicate consecutive vertex at edge $i"))
        for j in i+1:n
            # Adjacent edges share a legitimate endpoint.
            (j == i || j == mod1(i+1,n) || i == mod1(j+1,n)) && continue
            c, d = p[j], p[mod1(j+1,n)]
            _segments_intersect(a,b,c,d) &&
                throw(ArgumentError("polygon self-intersects at edges $i and $j"))
        end
    end
    return p
end
