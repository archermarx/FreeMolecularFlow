_temperature(b::Inflow) = b.temperature
_temperature(b::BackPressure) = b.temperature
_temperature(b::DiffuseWall) = b.temperature

function _prescribed_fluxes(mesh::RZMesh, gas::Gas)
    # q_j is the outward emission number flux from boundary segment j
    # [particles m⁻² s⁻¹]. Diffuse-wall entries remain zero here because their
    # emission is the unknown solved by `_solve_radiosity`.
    q = zeros(length(mesh.boundary_segments))
    label_areas = Dict{String,Float64}()
    for s in mesh.boundary_segments
        label_areas[s.label] = get(label_areas,s.label,0.0) + s.area
    end
    for (i,s) in pairs(mesh.boundary_segments)
        if s.condition isa Inflow
            # A single inflow specification is a total mass rate for its label.
            # Dividing by molecular mass and total labeled area makes the
            # emitted number flux invariant to boundary mesh refinement.
            q[i] = s.condition.mass_flow_rate /
                   (gas.molecular_mass * label_areas[s.label])
        elseif s.condition isa BackPressure
            # An isotropic Maxwellian strikes a plane with one-sided number
            # flux n*c̄/4. The same expression gives injection from a reservoir.
            n = s.condition.pressure / (BOLTZMANN*s.condition.temperature)
            q[i] = n * mean_molecular_speed(gas,s.condition.temperature) / 4
        end
    end
    q
end

function _solve_radiosity(mesh::RZMesh, H, gas::Gas, tolerance)
    segments = mesh.boundary_segments
    diffuse = findall(s -> s.condition isa DiffuseWall, segments)
    prescribed = findall(s -> !(s.condition isa DiffuseWall), segments)
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

function _evaluate_fields(mesh::RZMesh, sample_points, patches, bvh, q, gas, tol;
                          reporter=nothing, closure=NaN, radiosity=NaN,
                          phase=:field_reconstruction)
    labels = sort(unique(s.label for s in mesh.boundary_segments))
    count = length(sample_points)
    direct = Dict(label => zeros(count) for label in labels)
    contribution = Dict(label => zeros(count) for label in labels)
    momentum = zeros(3,count)
    density = zeros(count)
    speeds = [mean_molecular_speed(gas,_temperature(s.condition)) for s in mesh.boundary_segments]

    # Axisymmetry permits all field points to be placed at Cartesian azimuth 0.
    # The full revolved source geometry still contributes from every azimuth.
    for (ic,rz) in pairs(sample_points)
        point = (rz[1],rz[2],0.0)
        for (ip,patch) in pairs(patches)
            ray = _vsub(point,patch.center)
            distance = _norm(ray)
            distance > tol || continue
            direction = _vscale(ray,inv(distance))
            _dot(patch.normal,direction) > 0 || continue
            _occluded(bvh,patch.center,point,ip,0,tol) && continue
            omega = _solid_angle(point,patch) * _area_scale(patch)
            omega > 0 || continue
            iseg = patch.segment
            label = mesh.boundary_segments[iseg].label
            # Katz's direct view factor is the visible solid-angle fraction of
            # the full sphere. It is geometric and independent of source flux.
            direct[label][ic] += omega/(4pi)
            # For cosine-law Maxwellian emission, q = n_source*c̄/4 and the
            # local density contribution integrates to q*Ω/(π*c̄).
            dn = q[iseg] * omega / (pi*speeds[iseg])
            contribution[label][ic] += dn
            density[ic] += dn
            # Integral of molecular direction over the spherical triangle.
            # `_vector_solid_angle` points from the field point to the source,
            # hence the minus sign for source-to-field molecular motion.
            vector_omega = _vscale(_vector_solid_angle(point,patch),-_area_scale(patch))
            momentum[1,ic] += q[iseg]/pi * vector_omega[1]
            momentum[2,ic] += q[iseg]/pi * vector_omega[2]
            momentum[3,ic] += q[iseg]/pi * vector_omega[3]
        end
        if density[ic] > 100eps(Float64)
            momentum[:,ic] ./= density[ic]
        else
            momentum[:,ic] .= 0
        end
        momentum[3,ic] = 0.0 # exact by axisymmetry; removes finite-wedge roundoff
        _status!(reporter,phase,ic;total=count,
                 exchange_closure=closure,radiosity=radiosity)
    end
    _status!(reporter,phase,count;total=count,exchange_closure=closure,
             radiosity=radiosity,force=true)
    labels,direct,contribution,density,momentum
end

function _cell_fields(mesh::RZMesh, patches, bvh, q, gas, tol,
                      reporter=nothing, closure=NaN, radiosity=NaN)
    _evaluate_fields(mesh,mesh.centers,patches,bvh,q,gas,tol;
                     reporter,closure,radiosity)
end

function _particle_balance(mesh,H,q)
    # Diffuse walls have zero net particle transfer by construction. Therefore
    # steady global balance compares only prescribed emission through openings
    # with the incident flux that ultimately escapes through those openings.
    openings = findall(s -> !(s.condition isa DiffuseWall), mesh.boundary_segments)
    prescribed_emission = sum(mesh.boundary_segments[i].area*q[i] for i in openings)
    incident = H*q
    escaped = sum(mesh.boundary_segments[i].area*incident[i] for i in openings)
    abs(prescribed_emission-escaped)/max(prescribed_emission,escaped,1.0)
end

"""
    solve(geometry, boundaries, gas; options=SolverOptions(),
          status_interval=0, status_io=stdout) -> FlowResult

Solve steady, single-species, collisionless free-molecular flow in an
axisymmetric domain. Boundary labels in `geometry` map to `Inflow`,
`BackPressure`, `DiffuseWall`, or `Axis` values in `boundaries`.
"""
function solve(geometry::AxisymmetricGeometry, boundaries, gas::Gas;
               options::SolverOptions=SolverOptions(),
               status_interval::Real=0.0, status_io::IO=stdout,
               status_reporter::Union{Nothing,StatusReporter}=nothing)
    reporter = status_reporter === nothing && status_interval > 0 ?
               StatusReporter(status_interval,status_io) : status_reporter
    status_interval >= 0 || throw(ArgumentError("status_interval must be nonnegative"))
    _status!(reporter,:mesh,0;force=true)
    # 1. Build the R-Z solution mesh and its refined physical boundary rings.
    mesh = _make_mesh(geometry,boundaries,options)
    _status!(reporter,:mesh,length(mesh.cells);total=length(mesh.cells),force=true)
    _status!(reporter,:surface_quadrature,0;force=true)
    # 2. Revolve rings into a 3-D quadrature surface and acceleration structure.
    patches = _revolve_segments(mesh.boundary_segments,options.azimuthal_divisions)
    bvh = _build_bvh(patches)
    _status!(reporter,:surface_quadrature,length(patches);total=length(patches),force=true)
    scale = maximum(abs,Iterators.flatten(geometry.points);init=1.0)
    tol = options.visibility_tolerance * max(scale,1.0)
    # 3. Precompute the boundary-to-boundary transport operator H.
    H, closure = _boundary_exchange(mesh,patches,bvh,tol,
                                    options.azimuthal_divisions,reporter)
    _status!(reporter,:radiosity,0;exchange_closure=closure,force=true)
    # 4. Couple prescribed sources to conservative diffuse-wall re-emission.
    q, residual = _solve_radiosity(mesh,H,gas,options.radiosity_tolerance)
    _status!(reporter,:radiosity,1;total=1,exchange_closure=closure,
             radiosity=residual,force=true)
    # 5. Reconstruct geometric factors and gas moments at every cell center.
    labels,direct,contributions,density,velocity =
        _cell_fields(mesh,patches,bvh,q,gas,tol,reporter,closure,residual)
    balance = _particle_balance(mesh,H,q)
    _status!(reporter,:complete,1;total=1,exchange_closure=closure,
             radiosity=residual,particle_balance=balance,force=true)
    FlowResult(mesh,labels,direct,contributions,density,velocity,q,
               residual,balance,closure,gas,options)
end
