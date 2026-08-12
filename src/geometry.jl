function _normalise_boundaries(boundaries)
    Dict{String,BoundaryCondition}(String(k) => v for (k,v) in pairs(boundaries))
end

function _validate_options(options::SolverOptions)
    options.max_area >= 0 || throw(ArgumentError("max_area must be nonnegative (zero selects an automatic value)"))
    0 < options.min_angle < 34 || throw(ArgumentError("min_angle must lie in (0,34) degrees"))
    options.azimuthal_divisions >= 4 || throw(ArgumentError("azimuthal_divisions must be at least 4"))
    options.radiosity_tolerance > 0 || throw(ArgumentError("radiosity_tolerance must be positive"))
    options.visibility_tolerance > 0 || throw(ArgumentError("visibility_tolerance must be positive"))
    options.max_mesh_points >= 3 || throw(ArgumentError("max_mesh_points must be at least 3"))
end

function _original_edge_index(point_a, point_b, geometry::AxisymmetricGeometry)
    scale = maximum(abs, Iterators.flatten(geometry.points); init=1.0)
    atol = 2e-9 * max(scale, 1.0)
    for i in eachindex(geometry.points)
        a, b = geometry.points[i], geometry.points[mod1(i+1, length(geometry.points))]
        if _on_segment(a, b, point_a; atol) && _on_segment(a, b, point_b; atol)
            return i
        end
    end
    throw(ErrorException("mesh refinement created a boundary edge that cannot be mapped to the input polygon"))
end

function _triangle_area(a,b,c)
    abs(_orient(a,b,c)) / 2
end

function _make_mesh(geometry::AxisymmetricGeometry, boundaries, options::SolverOptions)
    _validate_options(options)
    bc = _normalise_boundaries(boundaries)
    used_labels = unique(geometry.edge_labels)
    missing = setdiff(used_labels, collect(keys(bc)))
    isempty(missing) || throw(ArgumentError("missing boundary conditions for labels: $(join(missing, ", "))"))
    extras = setdiff(collect(keys(bc)), used_labels)
    isempty(extras) || throw(ArgumentError("boundary conditions supplied for unused labels: $(join(extras, ", "))"))

    for i in eachindex(geometry.points)
        a, b = geometry.points[i], geometry.points[mod1(i+1,length(geometry.points))]
        condition = bc[geometry.edge_labels[i]]
        on_axis = a[2] == 0 && b[2] == 0
        if on_axis && !(condition isa Axis)
            throw(ArgumentError("edge $i lies on r=0 and must use Axis()"))
        elseif !on_axis && condition isa Axis
            throw(ArgumentError("Axis() boundary $(geometry.edge_labels[i]) must lie entirely on r=0"))
        end
    end

    points = copy(geometry.points)
    boundary_nodes = [collect(1:length(points)); 1]
    tri = triangulate(points; boundary_nodes, randomise=false, delete_ghosts=false)
    domain_area = abs(_signed_area(geometry.points))
    max_area = options.max_area > 0 ? options.max_area : domain_area / 200
    refine!(tri; min_angle=options.min_angle, max_area,
            max_points=options.max_mesh_points)

    mesh_points = NTuple{2,Float64}[(Float64(p[1]), Float64(p[2])) for p in get_points(tri)]
    cells = NTuple{3,Int}[]
    for T in each_solid_triangle(tri)
        t = (Int(T[1]), Int(T[2]), Int(T[3]))
        a,b,c = mesh_points[t[1]], mesh_points[t[2]], mesh_points[t[3]]
        _orient(a,b,c) < 0 && (t = (t[1],t[3],t[2]))
        push!(cells, t)
    end
    isempty(cells) && throw(ErrorException("mesher produced no interior cells"))

    centers = NTuple{2,Float64}[]
    volumes = Float64[]
    for t in cells
        a,b,c = mesh_points[t[1]], mesh_points[t[2]], mesh_points[t[3]]
        center = ((a[1]+b[1]+c[1])/3, (a[2]+b[2]+c[2])/3)
        push!(centers, center)
        push!(volumes, 2pi * center[2] * _triangle_area(a,b,c))
    end

    segments = BoundarySegment[]
    for e in each_segment(tri)
        a, b = mesh_points[Int(e[1])], mesh_points[Int(e[2])]
        original = _original_edge_index(a, b, geometry)
        oa, ob = geometry.points[original], geometry.points[mod1(original+1,length(geometry.points))]
        if (b[1]-a[1])*(ob[1]-oa[1]) + (b[2]-a[2])*(ob[2]-oa[2]) < 0
            a, b = b, a
        end
        label = geometry.edge_labels[original]
        condition = bc[label]
        slant = hypot(b[1]-a[1], b[2]-a[2])
        area = pi * (a[2] + b[2]) * slant
        area <= 10eps(Float64) && continue # symmetry axis
        push!(segments, BoundarySegment(a,b,label,condition,area))
    end
    isempty(segments) && throw(ArgumentError("geometry has no nonzero-area physical boundaries"))
    sort!(segments; by=s -> (s.label, s.a, s.b))
    return RZMesh(mesh_points, cells, centers, volumes, segments)
end
