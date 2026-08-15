"""
Geometry-dependent transport data that can be reused across gas and boundary
parameter sweeps. Construct with [`prepare`](@ref), then call `solve` one or
more times without rebuilding the mesh, surface quadrature, BVH, or exchange
matrix, or repeating cell-to-surface visibility queries.
"""
struct RadiositySystem
    diffuse_mask::BitVector
    diffuse::Vector{Int}
    prescribed::Vector{Int}
    matrix::Matrix{Float64}
    prescribed_coupling::Matrix{Float64}
    factorization::Union{Nothing,LU{Float64,Matrix{Float64},Vector{Int}}}
end

struct FieldAssemblyCache
    labels::Vector{String}
    segment_labels::Vector{Int}
    direct_view_factors::Dict{String,Vector{Float64}}
    label_areas::Dict{String,Float64}
end

struct PreparedSolver
    geometry::AxisymmetricGeometry
    mesh::RZMesh
    patches::Vector{SurfacePatch}
    occluder::AnalyticOccluder
    exchange_matrix::Matrix{Float64}
    exchange_closure_error::Float64
    solid_angles::Matrix{Float64}
    direction_moments::Array{Float64,3}
    second_direction_moments::Array{Float64,3}
    radiosity_system::RadiositySystem
    field_cache::FieldAssemblyCache
    azimuthal_divisions::Int
    azimuthal_convergence_error::Float64
    visibility_tolerance::Float64
    options::SolverOptions
end

struct DirectEvaluator
    patches::Vector{SurfacePatch}
    occluder::AnalyticOccluder
    tol::Float64
    azimuthal_divisions::Int
end

"""
Results of a steady free-molecular calculation.

All spatial arrays are cell-centered. Per-label dictionaries separate direct
geometric visibility from the actual density contribution after diffuse-wall
recycling. `temperature` contains the centered diagonal second moments in
kelvin, ordered `(T_z,T_r,T_θ)`; it does not imply a Maxwellian VDF. The saved
direct evaluator allows extraction lines to reuse the surface BVH without
retaining the larger parameter-sweep moment cache in every result.
"""
struct FlowResult
    mesh::RZMesh
    labels::Vector{String}
    direct_view_factors::Dict{String,Vector{Float64}}
    density_contributions::Dict{String,Vector{Float64}}
    density::Vector{Float64}
    velocity::Matrix{Float64}       # 3 × ncells, ordered (z,r,azimuthal)
    temperature::Matrix{Float64}    # 3 × ncells [K], centered (z,r,azimuthal)
    boundary_flux::Vector{Float64}  # particles m^-2 s^-1
    radiosity_residual::Float64
    particle_balance_residual::Float64
    exchange_closure_error::Float64
    azimuthal_divisions::Int
    azimuthal_convergence_error::Float64
    gas::Gas
    options::SolverOptions
    evaluator::DirectEvaluator
end

struct GeometricTransport
    patches::Vector{SurfacePatch}
    occluder::AnalyticOccluder
    exchange_matrix::Matrix{Float64}
    exchange_closure_error::Float64
    solid_angles::Matrix{Float64}
    direction_moments::Array{Float64,3}
    second_direction_moments::Array{Float64,3}
end

number_density(result::FlowResult) = result.density

_temperature(b::Inflow) = b.temperature
_temperature(b::BackPressure) = b.temperature
_temperature(b::DiffuseWall) = b.temperature

_is_diffuse(segment) = segment.condition isa DiffuseWall

function _prescribed_flux(boundary::Inflow,gas,label_area)
    # Inflow is a total mass rate shared by every segment with the same label.
    boundary.mass_flow_rate/(gas.molecular_mass*label_area)
end

function _prescribed_flux(boundary::BackPressure,gas,label_area)
    # A Maxwellian reservoir injects its one-sided impingement flux n*cbar/4.
    density = boundary.pressure/(BOLTZMANN*boundary.temperature)
    density*mean_molecular_speed(gas,boundary.temperature)/4
end

_prescribed_flux(::DiffuseWall,gas,label_area) = 0.0

function _label_areas(mesh::RZMesh)
    areas = Dict{String,Float64}()
    for segment in mesh.boundary_segments
        areas[segment.label] = get(areas,segment.label,0.0) + segment.area
    end
    areas
end

function _prescribed_fluxes(mesh::RZMesh,gas::Gas,label_areas=_label_areas(mesh))
    # q_j is the outward emission number flux from boundary segment j
    # [particles m⁻² s⁻¹]. Diffuse-wall entries remain zero here because their
    # emission is the unknown solved by `_solve_radiosity`.
    q = zeros(length(mesh.boundary_segments))
    for (i,segment) in pairs(mesh.boundary_segments)
        q[i] = _prescribed_flux(segment.condition,gas,label_areas[segment.label])
    end
    q
end

function _prepare_radiosity_system(mesh::RZMesh,H,tolerance)
    segments = mesh.boundary_segments
    mask = BitVector(_is_diffuse.(segments))
    diffuse, prescribed = findall(mask), findall(.!mask)
    isempty(prescribed) && return RadiositySystem(
        mask,diffuse,prescribed,zeros(0,0),zeros(0,0),nothing)
    isempty(diffuse) && return RadiositySystem(
        mask,diffuse,prescribed,zeros(0,0),zeros(0,length(prescribed)),nothing)
    # A diffuse wall emits exactly what it receives:
    #   q_D = H_DD q_D + H_DP q_P.
    # Rearranging gives one dense, deterministic linear radiosity solve.
    M = Matrix(I,length(diffuse),length(diffuse)) - H[diffuse,diffuse]
    # A perfectly closed diffuse enclosure has eigenvalue one in H_DD, so its
    # absolute density is arbitrary. Detect that physical non-uniqueness before
    # relying on a numerically unstable factorization.
    kappa = cond(M)
    (!isfinite(kappa) || kappa > inv(tolerance)) && throw(ArgumentError(
        "diffuse-wall radiosity system is singular or ill-conditioned; ensure the domain has a finite-area escape opening"))
    RadiositySystem(mask,diffuse,prescribed,M,H[diffuse,prescribed],lu(M))
end

function _solve_radiosity(mesh::RZMesh,gas::Gas,tolerance,
                          system::RadiositySystem,
                          label_areas=_label_areas(mesh))
    isempty(system.prescribed) && throw(ArgumentError(
        "a domain containing only diffuse walls has no unique steady density; add an inflow or back-pressure opening"))
    q = _prescribed_fluxes(mesh,gas,label_areas)
    isempty(system.diffuse) && return q, 0.0
    M = system.matrix
    rhs = system.prescribed_coupling * q[system.prescribed]
    qd = system.factorization \ rhs
    scale = max(maximum(abs,qd; init=0.0),maximum(abs,q; init=0.0),1.0)
    # Roundoff may create tiny negative values, but a materially negative
    # outward wall flux violates particle conservation and indicates failure.
    minimum(qd; init=0.0) < -100eps(Float64)*scale &&
        throw(ErrorException("radiosity solve produced a negative wall flux"))
    q[system.diffuse] = max.(qd,0.0)
    residual = norm(M*q[system.diffuse]-rhs) / max(norm(rhs),1.0)
    residual <= 10tolerance || throw(ErrorException(
        "diffuse-wall solve did not converge: relative residual $residual"))
    q, residual
end

function _solve_radiosity(mesh::RZMesh,H,gas::Gas,tolerance)
    system = _prepare_radiosity_system(mesh,H,tolerance)
    _solve_radiosity(mesh,gas,tolerance,system)
end

function _field_patch_quadrature(patches,ntheta)
    if iseven(ntheta)
        # `_revolve_segments` mirrors triangle diagonals for even quadratures,
        # so positive- and negative-y patches are exact pairs for a field point
        # in the meridional plane. Integrating one member with weight two halves
        # the dominant ray count without making a center-ray approximation.
        positive = findall(p -> p.center[3] > 0.0,patches)
        length(positive)*2 == length(patches) && return positive,2.0
    end
    eachindex(patches),1.0
end

function _geometric_field_moments(mesh::RZMesh,sample_points,patches,bvh,tol;
                                  ntheta::Int=1,reporter=nothing,
                                  closure=NaN,radiosity=NaN,
                                  phase=:field_reconstruction)
    count = length(sample_points)
    segment_count = length(mesh.boundary_segments)
    solid_angles = zeros(segment_count,count)
    # Axisymmetry makes the azimuthal moment exactly zero. Cache only axial and
    # radial components, reducing persistent sweep storage by one quarter.
    direction_moments = zeros(2,segment_count,count)
    # Diagonal of ∫ ŝ ŝᵀ dΩ, ordered (z,r,azimuthal). Axisymmetry eliminates
    # the azimuthal cross moments; the diagonal determines directional
    # temperatures without assuming that the local VDF is Maxwellian.
    second_direction_moments = zeros(3,segment_count,count)
    patch_indices,quadrature_weight = _field_patch_quadrature(patches,ntheta)
    # Each worker owns a traversal stack. Reusing it makes `_occluded`
    # allocation-free while avoiding synchronization inside BVH traversal.
    traversal_stacks = [_new_traversal_stack(bvh) for _ in 1:Threads.maxthreadid()]
    completed = Threads.Atomic{Int}(0)
    reporter_lock = ReentrantLock()
    tolerance_squared = tol^2

    # Axisymmetry permits all field points to be placed at Cartesian azimuth 0.
    # Cells are independent, so the outer loop scales across Julia threads.
    Threads.@threads for ic in eachindex(sample_points)
        rz = sample_points[ic]
        point = (rz[1],rz[2],0.0)
        traversal_stack = traversal_stacks[Threads.threadid()]
        for ip in patch_indices
            patch = patches[ip]
            ray = _vsub(point,patch.center)
            _dot(ray,ray) > tolerance_squared || continue
            # Normalization cannot change the sign of the facing test.
            _dot(patch.normal,ray) > 0 || continue
            _occluded(bvh,patch.surface_center,point,ip,0,tol,
                      traversal_stack) && continue
            omega,vector_omega,second_omega = _solid_angle_moments(point,patch)
            scale = quadrature_weight*patch.area_scale
            omega *= scale
            omega > 0 || continue
            iseg = patch.segment
            solid_angles[iseg,ic] += omega
            # Integral of molecular direction over the spherical triangle.
            # The vector moment points from the field point to the source,
            # hence the minus sign for source-to-field molecular motion.
            vector_omega = _vscale(vector_omega,-scale)
            direction_moments[1,iseg,ic] += vector_omega[1]
            direction_moments[2,iseg,ic] += vector_omega[2]
            for component in 1:3
                second_direction_moments[component,iseg,ic] +=
                    scale*second_omega[component]
            end
        end
        _thread_status!(reporter,reporter_lock,completed,phase,count;
                        exchange_closure=closure,radiosity=radiosity)
    end
    _status!(reporter,phase,count;total=count,exchange_closure=closure,
             radiosity=radiosity,force=true)
    solid_angles,direction_moments,second_direction_moments
end

function _field_assembly_cache(mesh::RZMesh,solid_angles)
    labels = sort(unique(s.label for s in mesh.boundary_segments))
    count = size(solid_angles,2)
    label_index = Dict(label=>i for (i,label) in pairs(labels))
    segment_label = [label_index[s.label] for s in mesh.boundary_segments]
    direct_values = zeros(length(labels),count)
    for ic in 1:count, iseg in eachindex(mesh.boundary_segments)
        direct_values[segment_label[iseg],ic] += solid_angles[iseg,ic]/(4pi)
    end
    direct = Dict(label=>collect(@view direct_values[i,:])
                  for (i,label) in pairs(labels))
    FieldAssemblyCache(labels,segment_label,direct,_label_areas(mesh))
end

function _assemble_fields(mesh::RZMesh,solid_angles,direction_moments,
                          second_direction_moments,q,gas,
                          cache=_field_assembly_cache(mesh,solid_angles))
    labels = cache.labels
    count = size(solid_angles,2)
    contribution_values = zeros(length(labels),count)
    momentum = zeros(3,count)
    raw_second_moment = zeros(3,count)
    temperature = zeros(3,count)
    density = zeros(count)
    speeds = [mean_molecular_speed(gas,_temperature(s.condition))
              for s in mesh.boundary_segments]
    density_weights = q ./ (pi .* speeds)
    momentum_weights = q ./ pi
    second_moment_weights = density_weights .* [
        3BOLTZMANN*_temperature(s.condition)/gas.molecular_mass
        for s in mesh.boundary_segments]
    for ic in 1:count, iseg in eachindex(mesh.boundary_segments)
        omega = solid_angles[iseg,ic]
        omega > 0 || continue
        ilabel = cache.segment_labels[iseg]
        # For cosine-law Maxwellian emission, q = n_source*c̄/4 and the local
        # density contribution is q*Ω/(π*c̄).
        dn = density_weights[iseg]*omega
        contribution_values[ilabel,ic] += dn
        density[ic] += dn
        momentum[1,ic] += momentum_weights[iseg]*direction_moments[1,iseg,ic]
        momentum[2,ic] += momentum_weights[iseg]*direction_moments[2,iseg,ic]
        for component in 1:3
            raw_second_moment[component,ic] += second_moment_weights[iseg] *
                second_direction_moments[component,iseg,ic]
        end
    end
    for ic in 1:count
        if density[ic] > 100eps(Float64)
            momentum[1,ic] /= density[ic]
            momentum[2,ic] /= density[ic]
            for component in 1:3
                variance = raw_second_moment[component,ic]/density[ic] -
                           momentum[component,ic]^2
                # The analytic moments make the covariance nonnegative. Clamp
                # only the possible negative zero from floating-point roundoff.
                temperature[component,ic] =
                    gas.molecular_mass/BOLTZMANN * max(variance,0.0)
            end
        else
            momentum[1,ic] = 0.0
            momentum[2,ic] = 0.0
        end
    end
    # Convert dense per-label contributions to the public dictionary form only
    # after the numerical loop.
    contribution = Dict(label=>collect(@view contribution_values[i,:])
                        for (i,label) in pairs(labels))
    # Results own their mutable arrays; callers cannot corrupt PreparedSolver's
    # cached geometry by modifying a previous FlowResult.
    direct = Dict(label=>copy(values)
                  for (label,values) in cache.direct_view_factors)
    copy(labels),direct,contribution,density,momentum,temperature
end

function _evaluate_fields(mesh::RZMesh,sample_points,patches,bvh,q,gas,tol;
                          ntheta::Int=1,reporter=nothing,closure=NaN,
                          radiosity=NaN,phase=:field_reconstruction)
    solid_angles,direction_moments,second_direction_moments = _geometric_field_moments(
        mesh,sample_points,patches,bvh,tol;
        ntheta,reporter,closure,radiosity,phase)
    _assemble_fields(mesh,solid_angles,direction_moments,
                     second_direction_moments,q,gas)
end

function _particle_balance(mesh,H,q)
    # Diffuse walls have zero net particle transfer by construction. Therefore
    # steady global balance compares only prescribed emission through openings
    # with the incident flux that ultimately escapes through those openings.
    openings = findall(s -> !_is_diffuse(s),mesh.boundary_segments)
    prescribed_emission = sum(mesh.boundary_segments[i].area*q[i] for i in openings)
    incident = H*q
    escaped = sum(mesh.boundary_segments[i].area*incident[i] for i in openings)
    abs(prescribed_emission-escaped)/max(prescribed_emission,escaped,1.0)
end

function _make_reporter(status_interval,status_io,status_reporter)
    status_interval >= 0 || throw(ArgumentError("status_interval must be nonnegative"))
    status_reporter === nothing && status_interval > 0 ?
        StatusReporter(status_interval,status_io) : status_reporter
end

function _geometric_transport(geometry,mesh,nazimuth,tol,reporter)
    _status!(reporter,:surface_quadrature,0;force=true)
    patches = _revolve_segments(mesh.boundary_segments,nazimuth)
    occluder = _build_occluder(geometry)
    _status!(reporter,:surface_quadrature,length(patches);
             total=length(patches),force=true)
    H,closure = _boundary_exchange(
        mesh,patches,occluder,tol,nazimuth,reporter)
    solid_angles,direction_moments,second_direction_moments = _geometric_field_moments(
        mesh,mesh.centers,patches,occluder,tol;
        ntheta=nazimuth,reporter,closure,phase=:field_precompute)
    GeometricTransport(patches,occluder,H,closure,solid_angles,direction_moments,
                       second_direction_moments)
end

function _relative_change(current,previous)
    difference_scale = 0.0
    value_scale = eps(Float64)
    for i in eachindex(current,previous)
        difference_scale = max(difference_scale,abs(current[i]-previous[i]))
        value_scale = max(value_scale,abs(current[i]),abs(previous[i]))
    end
    difference_scale/value_scale
end

function _transport_change(current::GeometricTransport,previous::GeometricTransport)
    max(_relative_change(current.exchange_matrix,previous.exchange_matrix),
        _relative_change(current.solid_angles,previous.solid_angles),
        _relative_change(current.direction_moments,previous.direction_moments),
        _relative_change(current.second_direction_moments,
                         previous.second_direction_moments))
end

function _converged_transport(geometry,mesh,options,tol,reporter)
    divisions = options.azimuthal_divisions
    transport = _geometric_transport(geometry,mesh,divisions,tol,reporter)
    options.azimuthal_tolerance == 0 && return transport,divisions,NaN

    error = Inf
    while divisions < options.max_azimuthal_divisions
        next_divisions = min(2divisions,options.max_azimuthal_divisions)
        candidate = _geometric_transport(
            geometry,mesh,next_divisions,tol,reporter)
        error = _transport_change(candidate,transport)
        _status!(reporter,:azimuthal_convergence,next_divisions;
                 total=options.max_azimuthal_divisions,
                 exchange_closure=error,force=true)
        error <= options.azimuthal_tolerance &&
            return candidate,next_divisions,error
        transport,divisions = candidate,next_divisions
    end
    throw(ErrorException(
        "azimuthal quadrature did not reach tolerance $(options.azimuthal_tolerance) " *
        "by max_azimuthal_divisions=$(options.max_azimuthal_divisions) " *
        "(last relative change: $error)"))
end

"""
    prepare(geometry, boundaries; options=SolverOptions(),
            status_interval=0, status_io=stdout) -> PreparedSolver

Construct all geometry-dependent transport data. The returned object can be
reused when only gas properties or boundary conditions change.
"""
function prepare(geometry::AxisymmetricGeometry,boundaries;
                 options::SolverOptions=SolverOptions(),
                 status_interval::Real=0.0,status_io::IO=stdout,
                 status_reporter::Union{Nothing,StatusReporter}=nothing)
    reporter = _make_reporter(status_interval,status_io,status_reporter)
    _status!(reporter,:mesh,0;force=true)
    mesh = _make_mesh(geometry,boundaries,options)
    _status!(reporter,:mesh,length(mesh.cells);total=length(mesh.cells),force=true)
    scale = maximum(abs,Iterators.flatten(geometry.points);init=1.0)
    tol = options.visibility_tolerance*max(scale,1.0)
    transport,nazimuth,azimuthal_error = _converged_transport(
        geometry,mesh,options,tol,reporter)
    H,closure = transport.exchange_matrix,transport.exchange_closure_error
    solid_angles = transport.solid_angles
    direction_moments = transport.direction_moments
    second_direction_moments = transport.second_direction_moments
    radiosity_system = _prepare_radiosity_system(
        mesh,H,options.radiosity_tolerance)
    field_cache = _field_assembly_cache(mesh,solid_angles)
    PreparedSolver(geometry,mesh,transport.patches,transport.occluder,H,closure,solid_angles,
                   direction_moments,second_direction_moments,
                   radiosity_system,field_cache,nazimuth,
                   azimuthal_error,tol,options)
end

function _mesh_with_boundaries(prepared::PreparedSolver,boundaries)
    bc = _validate_boundaries(prepared.geometry,boundaries)
    old = prepared.mesh
    segments = [BoundarySegment(s.a,s.b,s.label,bc[s.label],s.area)
                for s in old.boundary_segments]
    RZMesh(old.points,old.cells,old.centers,old.volumes,segments)
end

function _solve_prepared(prepared::PreparedSolver,mesh::RZMesh,gas::Gas,reporter)
    closure = prepared.exchange_closure_error
    H = prepared.exchange_matrix
    options = prepared.options
    _status!(reporter,:radiosity,0;exchange_closure=closure,force=true)
    cached_system = prepared.radiosity_system
    same_pattern = all(i -> cached_system.diffuse_mask[i] ==
                            _is_diffuse(mesh.boundary_segments[i]),
                       eachindex(mesh.boundary_segments))
    system = same_pattern ? cached_system :
             _prepare_radiosity_system(mesh,H,options.radiosity_tolerance)
    q,residual = _solve_radiosity(
        mesh,gas,options.radiosity_tolerance,system,
        prepared.field_cache.label_areas)
    _status!(reporter,:radiosity,1;total=1,exchange_closure=closure,
             radiosity=residual,force=true)
    labels,direct,contributions,density,velocity,temperature = _assemble_fields(
        mesh,prepared.solid_angles,prepared.direction_moments,
        prepared.second_direction_moments,q,gas,
        prepared.field_cache)
    balance = _particle_balance(mesh,H,q)
    _status!(reporter,:complete,1;total=1,exchange_closure=closure,
             radiosity=residual,particle_balance=balance,force=true)
    FlowResult(mesh,labels,direct,contributions,density,velocity,temperature,q,
               residual,balance,closure,prepared.azimuthal_divisions,
               prepared.azimuthal_convergence_error,gas,options,
               DirectEvaluator(prepared.patches,prepared.occluder,
                               prepared.visibility_tolerance,
                               prepared.azimuthal_divisions))
end

"""Solve using the boundary conditions stored when `prepared` was constructed."""
function solve(prepared::PreparedSolver,gas::Gas;
               status_interval::Real=0.0,status_io::IO=stdout,
               status_reporter::Union{Nothing,StatusReporter}=nothing)
    reporter = _make_reporter(status_interval,status_io,status_reporter)
    _solve_prepared(prepared,prepared.mesh,gas,reporter)
end

"""
Solve with replacement boundary conditions while reusing a prepared geometry.
Labels and geometry must match, but boundary values and boundary types may
change.
"""
function solve(prepared::PreparedSolver,boundaries,gas::Gas;
               status_interval::Real=0.0,status_io::IO=stdout,
               status_reporter::Union{Nothing,StatusReporter}=nothing)
    reporter = _make_reporter(status_interval,status_io,status_reporter)
    mesh = _mesh_with_boundaries(prepared,boundaries)
    _solve_prepared(prepared,mesh,gas,reporter)
end

"""
    solve(geometry, boundaries, gas; options=SolverOptions(),
          status_interval=0, status_io=stdout) -> FlowResult

Solve steady, single-species, collisionless free-molecular flow in an
axisymmetric domain. This convenience method prepares the transport geometry
and immediately solves it; use `prepare` explicitly for parameter sweeps.
"""
function solve(geometry::AxisymmetricGeometry,boundaries,gas::Gas;
               options::SolverOptions=SolverOptions(),
               status_interval::Real=0.0,status_io::IO=stdout,
               status_reporter::Union{Nothing,StatusReporter}=nothing)
    reporter = _make_reporter(status_interval,status_io,status_reporter)
    prepared = prepare(geometry,boundaries;options,status_reporter=reporter)
    solve(prepared,gas;status_reporter=reporter)
end
