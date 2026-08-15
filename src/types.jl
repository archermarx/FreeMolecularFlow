# Exact SI value of k_B and the CODATA atomic-mass conversion used by Gas.
const BOLTZMANN = 1.380649e-23
const ATOMIC_MASS = 1.66053906660e-27

abstract type BoundaryCondition end

"""
A prescribed cosine-law inflow with total mass flow rate (kg/s).

All refined boundary segments sharing the same label divide this total flow in
proportion to their revolved area; `mass_flow_rate` is therefore not a local
flux density.
"""
struct Inflow <: BoundaryCondition
    mass_flow_rate::Float64
    temperature::Float64
    function Inflow(mass_flow_rate::Real, temperature::Real)
        mass_flow_rate >= 0 || throw(ArgumentError("inflow mass_flow_rate must be nonnegative"))
        temperature > 0 || throw(ArgumentError("inflow temperature must be positive"))
        new(Float64(mass_flow_rate), Float64(temperature))
    end
end

"""
An open boundary connected to a stationary reservoir at pressure (Pa) and
temperature (K). Molecules may leave freely, while the reservoir supplies the
incoming half of an isotropic Maxwellian. Set pressure to zero for vacuum.
"""
struct BackPressure <: BoundaryCondition
    pressure::Float64
    temperature::Float64
    function BackPressure(pressure::Real, temperature::Real)
        pressure >= 0 || throw(ArgumentError("back-pressure pressure must be nonnegative"))
        temperature > 0 || throw(ArgumentError("back-pressure temperature must be positive"))
        new(Float64(pressure), Float64(temperature))
    end
end

"""
A fully accommodating, cosine-law diffuse wall.

Incident molecules lose memory of their incoming velocity and are re-emitted
at the wall temperature. The radiosity solve enforces zero net particle loss at
this boundary.
"""
struct DiffuseWall <: BoundaryCondition
    temperature::Float64
    function DiffuseWall(temperature::Real)
        temperature > 0 || throw(ArgumentError("wall temperature must be positive"))
        new(Float64(temperature))
    end
end

"""
A zero-area symmetry-axis edge. Axis edges must lie at `r = 0`; they close the
R-Z polygon but do not create a physical surface when revolved.
"""
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
    max_area::Float64 = 0.0             # Maximum R-Z triangle area [m²]; 0 is automatic.
    max_boundary_length::Float64 = 0.0  # Maximum transport-segment length [m]; 0 is automatic.
    min_angle::Float64 = 20.0           # Delaunay mesh-quality target [degrees].
    azimuthal_divisions::Int = 64       # Fixed resolution, or adaptive starting resolution.
    azimuthal_tolerance::Float64 = 0.0  # Relative convergence target; 0 disables adaptivity.
    max_azimuthal_divisions::Int = 256  # Adaptive-resolution cost guard.
    radiosity_tolerance::Float64 = 1e-10 # Linear-solve conditioning/residual tolerance.
    visibility_tolerance::Float64 = 1e-10 # Relative ray-intersection tolerance.
    max_mesh_points::Int = 200_000      # Guard against runaway mesh refinement.
end

const EXTRACTION_FIELDS = (:number_density, :velocity, :temperature, :view_factors,
                           :density_contributions)

function _extraction_choice(value, name, allowed)
    choice = Symbol(lowercase(String(value)))
    choice in allowed || throw(ArgumentError(
        "extraction-line $name must be " * join("`" .* string.(allowed) .* "`", ", ")))
    choice
end

"""
A uniformly sampled straight or piecewise-linear path in the R-Z plane.

Exactly one of `num_points` and `spacing` must be supplied. `method` is
`:direct` or `:cell`; `outside_domain` is `:keep`, `:drop`, or `:error`.
"""
struct ExtractionLine
    name::String
    points::Vector{NTuple{2,Float64}}
    num_points::Union{Nothing,Int}
    spacing::Union{Nothing,Float64}
    method::Symbol
    outside_domain::Symbol
    fields::Vector{Symbol}
    filename::String
    function ExtractionLine(name::AbstractString, points;
                            num_points::Union{Nothing,Integer}=nothing,
                            spacing::Union{Nothing,Real}=nothing,
                            method::Union{Symbol,AbstractString}=:direct,
                            outside_domain::Union{Symbol,AbstractString}=:keep,
                            fields=EXTRACTION_FIELDS,
                            filename::Union{Nothing,AbstractString}=nothing)
        clean_name = strip(String(name))
        isempty(clean_name) && throw(ArgumentError("extraction-line name may not be empty"))
        pts = NTuple{2,Float64}[(Float64(p[1]),Float64(p[2])) for p in points]
        length(pts) >= 2 || throw(ArgumentError("extraction line needs at least two control points"))
        all(p -> all(isfinite,p),pts) ||
            throw(ArgumentError("extraction-line coordinates must be finite"))
        all(p -> p[2] >= 0,pts) ||
            throw(ArgumentError("extraction-line radius r must be nonnegative"))
        all(i -> pts[i] != pts[i+1],1:length(pts)-1) ||
            throw(ArgumentError("adjacent extraction-line control points must differ"))

        (isnothing(num_points) != isnothing(spacing)) || throw(ArgumentError(
            "extraction line requires exactly one of num_points or spacing"))
        count = isnothing(num_points) ? nothing : Int(num_points)
        step = isnothing(spacing) ? nothing : Float64(spacing)
        isnothing(count) || count >= 2 || throw(ArgumentError(
            "extraction-line num_points must be at least 2"))
        isnothing(step) || (isfinite(step) && step > 0) || throw(ArgumentError(
            "extraction-line spacing must be finite and positive"))

        method_symbol = _extraction_choice(method,"method",(:direct,:cell))
        outside_symbol = _extraction_choice(
            outside_domain,"outside_domain",(:keep,:drop,:error))
        selected = Symbol[Symbol(lowercase(String(field))) for field in fields]
        isempty(selected) && throw(ArgumentError("extraction-line fields may not be empty"))
        length(unique(selected)) == length(selected) || throw(ArgumentError(
            "extraction-line fields may not contain duplicates"))
        unknown = setdiff(selected,EXTRACTION_FIELDS)
        isempty(unknown) || throw(ArgumentError(
            "unknown extraction fields: $(join(string.(unknown), ", "))"))

        stem = strip(replace(lowercase(clean_name),r"[^a-z0-9_-]+" => "_"),'_')
        default_filename = (isempty(stem) ? "extraction" : stem) * ".csv"
        output = isnothing(filename) ? default_filename : String(filename)
        isempty(strip(output)) && throw(ArgumentError("extraction-line filename may not be empty"))
        new(clean_name,pts,count,step,method_symbol,outside_symbol,selected,output)
    end
end

mutable struct StatusReporter
    interval::Float64
    io::IO
    started::Float64
    last_printed::Float64
    header_printed::Bool
end

function StatusReporter(interval::Real, io::IO)
    interval >= 0 || throw(ArgumentError("status_interval must be nonnegative"))
    now = time()
    # Backdate last_printed so the initial phase row is emitted immediately.
    StatusReporter(Float64(interval),io,now,now-Float64(interval),false)
end

_residual_text(value) = isfinite(value) ? @sprintf("%.3e",value) : "-"

const _STATUS_PHASE_NAMES = Dict(
    :mesh => "mesh",
    :surface_quadrature => "surface",
    :boundary_exchange => "exchange",
    :exchange_balance => "balance",
    :radiosity => "radiosity",
    :field_precompute => "moments",
    :azimuthal_convergence => "azimuth",
    :field_reconstruction => "fields",
    :line_extraction => "extract",
    :complete => "complete")

function _status!(reporter::Union{Nothing,StatusReporter}, phase::Symbol,
                  iteration::Integer; total::Integer=0,
                  exchange_closure::Real=NaN, radiosity::Real=NaN,
                  particle_balance::Real=NaN, force::Bool=false)
    reporter === nothing && return
    now = time()
    # `force` is used at phase boundaries so even short calculations document
    # their progression; loop updates remain rate-limited by wall-clock time.
    (force || now-reporter.last_printed >= reporter.interval) || return
    progress = total > 0 ? "$(iteration)/$(total)" : string(iteration)
    if !reporter.header_printed
        @printf(reporter.io,"%-10s %11s %10s %11s %11s %11s\n",
                "phase","iteration","elapsed(s)","closure","radiosity","balance")
        @printf(reporter.io,"%-10s %11s %10s %11s %11s %11s\n",
                "----------","-----------","----------","-----------","-----------","-----------")
        reporter.header_printed = true
    end
    @printf(reporter.io,
        "%-10s %11s %10.2f %11s %11s %11s\n",
        get(_STATUS_PHASE_NAMES,phase,String(phase)),progress,now-reporter.started,
        _residual_text(exchange_closure),_residual_text(radiosity),
        _residual_text(particle_balance))
    flush(reporter.io)
    reporter.last_printed = now
end

function _thread_status!(reporter,mutex,completed,phase,total;kwargs...)
    reporter === nothing && return
    iteration = Threads.atomic_add!(completed,1) + 1
    lock(mutex) do
        _status!(reporter,phase,iteration;total,kwargs...)
    end
end

struct BoundarySegment
    # Endpoints follow the normalized counter-clockwise polygon direction. This
    # orientation makes the gas-side normal unambiguous during revolution.
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
    # Axisymmetric cell volume obtained by revolving each triangle around r=0.
    volumes::Vector{Float64}
    boundary_segments::Vector{BoundarySegment}
end

function _signed_area(p)
    # Shoelace sign establishes polygon orientation in the (z,r) plane.
    0.5 * sum(p[i][1] * p[mod1(i+1, length(p))][2] -
              p[mod1(i+1, length(p))][1] * p[i][2] for i in eachindex(p))
end

function _orient(a, b, c)
    # Twice the signed area of triangle (a,b,c); positive means a left turn.
    (b[1]-a[1]) * (c[2]-a[2]) - (b[2]-a[2]) * (c[1]-a[1])
end

function _on_segment(a, b, p; atol=1e-13)
    abs(_orient(a,b,p)) <= atol * max(1.0, hypot(b[1]-a[1], b[2]-a[2])) &&
        min(a[1],b[1])-atol <= p[1] <= max(a[1],b[1])+atol &&
        min(a[2],b[2])-atol <= p[2] <= max(a[2],b[2])+atol
end

function _segments_intersect(a,b,c,d)
    # Proper crossings are detected by opposing orientations; the final checks
    # include touching/collinear intersections, which are invalid for this
    # solver's single simple polygon.
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
