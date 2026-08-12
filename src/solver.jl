"""
Geometry-dependent transport data that can be reused across gas and boundary
parameter sweeps. Construct with [`prepare`](@ref), then call `solve` one or
more times without rebuilding the mesh, surface quadrature, BVH, or exchange
matrix, or repeating cell-to-surface visibility queries.
"""
struct PreparedSolver
    geometry::AxisymmetricGeometry
    mesh::RZMesh
    bvh::PatchBVH
    exchange_matrix::Matrix{Float64}
    exchange_closure_error::Float64
    solid_angles::Matrix{Float64}
    direction_moments::Array{Float64,3}
    visibility_tolerance::Float64
    options::SolverOptions
end

struct DirectEvaluator
    bvh::PatchBVH
    tol::Float64
end

"""
Results of a steady free-molecular calculation.

All spatial arrays are cell-centered. Per-label dictionaries separate direct
geometric visibility from the actual density contribution after diffuse-wall
recycling. The saved direct evaluator allows extraction lines to reuse the
surface BVH without retaining the larger parameter-sweep moment cache in every
result.
"""
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
    gas::Gas
    options::SolverOptions
    evaluator::DirectEvaluator
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

function _prescribed_fluxes(mesh::RZMesh, gas::Gas)
    # q_j is the outward emission number flux from boundary segment j
    # [particles m⁻² s⁻¹]. Diffuse-wall entries remain zero here because their
    # emission is the unknown solved by `_solve_radiosity`.
    q = zeros(length(mesh.boundary_segments))
    label_areas = Dict{String,Float64}()
    for s in mesh.boundary_segments
        label_areas[s.label] = get(label_areas,s.label,0.0) + s.area
    end
    for (i,segment) in pairs(mesh.boundary_segments)
        q[i] = _prescribed_flux(segment.condition,gas,label_areas[segment.label])
    end
    q
end

function _solve_radiosity(mesh::RZMesh, H, gas::Gas, tolerance)
    segments = mesh.boundary_segments
    diffuse = findall(_is_diffuse,segments)
    prescribed = findall(s -> !_is_diffuse(s),segments)
    isempty(prescribed) && throw(ArgumentError(
        "a domain containing only diffuse walls has no unique steady density; add an inflow or back-pressure opening"))
    q = _prescribed_fluxes(mesh,gas)
    isempty(diffuse) && return q, 0.0
    # A diffuse wall emits exactly what it receives:
    #   q_D = H_DD q_D + H_DP q_P.
    # Rearranging gives one dense, deterministic linear radiosity solve.
    M = Matrix(I,length(diffuse),length(diffuse)) - H[diffuse,diffuse]
    rhs = H[diffuse,prescribed] * q[prescribed]
    # A perfectly closed diffuse enclosure has eigenvalue one in H_DD, so its
    # absolute density is arbitrary. Detect that physical non-uniqueness before
    # relying on a numerically unstable factorization.
    kappa = cond(M)
    (!isfinite(kappa) || kappa > inv(tolerance)) && throw(ArgumentError(
        "diffuse-wall radiosity system is singular or ill-conditioned; ensure the domain has a finite-area escape opening"))
    qd = M \ rhs
    scale = max(maximum(abs,qd; init=0.0),maximum(abs,q; init=0.0),1.0)
    # Roundoff may create tiny negative values, but a materially negative
    # outward wall flux violates particle conservation and indicates failure.
    minimum(qd; init=0.0) < -100eps(Float64)*scale &&
        throw(ErrorException("radiosity solve produced a negative wall flux"))
    q[diffuse] = max.(qd,0.0)
    residual = norm(M*q[diffuse]-rhs) / max(norm(rhs),1.0)
    residual <= 10tolerance || throw(ErrorException(
        "diffuse-wall solve did not converge: relative residual $residual"))
    q, residual
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
            _occluded(bvh,patch.center,point,ip,0,tol,traversal_stack) && continue
            omega,vector_omega = _solid_angle_moments(point,patch)
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
        end
        _thread_status!(reporter,reporter_lock,completed,phase,count;
                        exchange_closure=closure,radiosity=radiosity)
    end
    _status!(reporter,phase,count;total=count,exchange_closure=closure,
             radiosity=radiosity,force=true)
    solid_angles,direction_moments
end

function _assemble_fields(mesh::RZMesh,solid_angles,direction_moments,q,gas)
    labels = sort(unique(s.label for s in mesh.boundary_segments))
    count = size(solid_angles,2)
    label_index = Dict(label=>i for (i,label) in pairs(labels))
    segment_label = [label_index[s.label] for s in mesh.boundary_segments]
    direct_values = zeros(length(labels),count)
    contribution_values = zeros(length(labels),count)
    momentum = zeros(3,count)
    density = zeros(count)
    speeds = [mean_molecular_speed(gas,_temperature(s.condition))
              for s in mesh.boundary_segments]
    density_weights = q ./ (pi .* speeds)
    momentum_weights = q ./ pi
    for ic in 1:count, iseg in eachindex(mesh.boundary_segments)
        omega = solid_angles[iseg,ic]
        omega > 0 || continue
        ilabel = segment_label[iseg]
        # Katz's direct view factor is the visible solid-angle fraction of the
        # full sphere and therefore needs no gas-dependent weighting.
        direct_values[ilabel,ic] += omega/(4pi)
        # For cosine-law Maxwellian emission, q = n_source*c̄/4 and the local
        # density contribution is q*Ω/(π*c̄).
        dn = density_weights[iseg]*omega
        contribution_values[ilabel,ic] += dn
        density[ic] += dn
        momentum[1,ic] += momentum_weights[iseg]*direction_moments[1,iseg,ic]
        momentum[2,ic] += momentum_weights[iseg]*direction_moments[2,iseg,ic]
    end
    for ic in 1:count
        if density[ic] > 100eps(Float64)
            momentum[1,ic] /= density[ic]
            momentum[2,ic] /= density[ic]
        else
            momentum[1,ic] = 0.0
            momentum[2,ic] = 0.0
        end
    end
    # Dictionaries remain the public label-oriented representation, but are
    # constructed only after the hot loop has finished.
    direct = Dict(label=>collect(@view direct_values[i,:])
                  for (i,label) in pairs(labels))
    contribution = Dict(label=>collect(@view contribution_values[i,:])
                        for (i,label) in pairs(labels))
    labels,direct,contribution,density,momentum
end

function _evaluate_fields(mesh::RZMesh,sample_points,patches,bvh,q,gas,tol;
                          ntheta::Int=1,reporter=nothing,closure=NaN,
                          radiosity=NaN,phase=:field_reconstruction)
    solid_angles,direction_moments = _geometric_field_moments(
        mesh,sample_points,patches,bvh,tol;
        ntheta,reporter,closure,radiosity,phase)
    _assemble_fields(mesh,solid_angles,direction_moments,q,gas)
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
    _status!(reporter,:surface_quadrature,0;force=true)
    patches = _revolve_segments(mesh.boundary_segments,options.azimuthal_divisions)
    bvh = _build_bvh(patches)
    _status!(reporter,:surface_quadrature,length(patches);total=length(patches),force=true)
    scale = maximum(abs,Iterators.flatten(geometry.points);init=1.0)
    tol = options.visibility_tolerance*max(scale,1.0)
    H,closure = _boundary_exchange(mesh,patches,bvh,tol,
                                   options.azimuthal_divisions,reporter)
    solid_angles,direction_moments = _geometric_field_moments(
        mesh,mesh.centers,patches,bvh,tol;
        ntheta=options.azimuthal_divisions,reporter,closure,
        phase=:field_precompute)
    PreparedSolver(geometry,mesh,bvh,H,closure,solid_angles,
                   direction_moments,tol,options)
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
    q,residual = _solve_radiosity(mesh,H,gas,options.radiosity_tolerance)
    _status!(reporter,:radiosity,1;total=1,exchange_closure=closure,
             radiosity=residual,force=true)
    labels,direct,contributions,density,velocity = _assemble_fields(
        mesh,prepared.solid_angles,prepared.direction_moments,q,gas)
    balance = _particle_balance(mesh,H,q)
    _status!(reporter,:complete,1;total=1,exchange_closure=closure,
             radiosity=residual,particle_balance=balance,force=true)
    FlowResult(mesh,labels,direct,contributions,density,velocity,q,
               residual,balance,closure,gas,options,
               DirectEvaluator(prepared.bvh,prepared.visibility_tolerance))
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
