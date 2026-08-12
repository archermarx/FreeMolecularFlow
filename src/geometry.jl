function _normalise_boundaries(boundaries)
    # Accept Dicts and NamedTuple-like mappings while keeping the numerical core
    # type-stable on String keys and the abstract boundary-condition interface.
    Dict{String,BoundaryCondition}(String(k) => v for (k,v) in pairs(boundaries))
end

function _validate_options(options::SolverOptions)
    isfinite(options.max_area) && options.max_area >= 0 ||
        throw(ArgumentError("max_area must be finite and nonnegative (zero selects an automatic value)"))
    isfinite(options.max_boundary_length) && options.max_boundary_length >= 0 ||
        throw(ArgumentError("max_boundary_length must be finite and nonnegative (zero selects an automatic value)"))
    0 < options.min_angle < 34 || throw(ArgumentError("min_angle must lie in (0,34) degrees"))
    options.azimuthal_divisions >= 4 || throw(ArgumentError("azimuthal_divisions must be at least 4"))
    isfinite(options.azimuthal_tolerance) && options.azimuthal_tolerance >= 0 ||
        throw(ArgumentError("azimuthal_tolerance must be finite and nonnegative"))
    options.max_azimuthal_divisions >= 4 ||
        throw(ArgumentError("max_azimuthal_divisions must be at least 4"))
    if options.azimuthal_tolerance > 0
        iseven(options.azimuthal_divisions) || throw(ArgumentError(
            "adaptive azimuthal quadrature requires an even azimuthal_divisions value"))
        iseven(options.max_azimuthal_divisions) || throw(ArgumentError(
            "adaptive azimuthal quadrature requires an even max_azimuthal_divisions value"))
        options.max_azimuthal_divisions >= options.azimuthal_divisions || throw(ArgumentError(
            "max_azimuthal_divisions must not be less than azimuthal_divisions"))
        options.max_azimuthal_divisions > options.azimuthal_divisions || throw(ArgumentError(
            "adaptive azimuthal quadrature requires max_azimuthal_divisions > azimuthal_divisions"))
    end
    options.radiosity_tolerance > 0 || throw(ArgumentError("radiosity_tolerance must be positive"))
    options.visibility_tolerance > 0 || throw(ArgumentError("visibility_tolerance must be positive"))
    options.max_mesh_points >= 3 || throw(ArgumentError("max_mesh_points must be at least 3"))
end

function _triangle_area(a,b,c)
    abs(_orient(a,b,c)) / 2
end

function _validate_boundaries(geometry::AxisymmetricGeometry,boundaries)
    bc = _normalise_boundaries(boundaries)
    labels = unique(geometry.edge_labels)
    supplied = collect(keys(bc))
    missing = setdiff(labels,supplied)
    isempty(missing) || throw(ArgumentError("missing boundary conditions for labels: $(join(missing, ", "))"))
    extras = setdiff(supplied,labels)
    isempty(extras) || throw(ArgumentError("boundary conditions supplied for unused labels: $(join(extras, ", "))"))

    # The axis is a topological boundary of the 2-D polygon but has zero area in
    # 3-D. Requiring Axis() explicitly prevents accidental sources at r=0.
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
    bc
end

function _make_boundary_segments(geometry::AxisymmetricGeometry,bc,
                                 max_boundary_length::Float64)
    segments = BoundarySegment[]
    for i in eachindex(geometry.points)
        parent_a = geometry.points[i]
        parent_b = geometry.points[mod1(i+1,length(geometry.points))]
        label = geometry.edge_labels[i]
        condition = bc[label]
        condition isa Axis && continue
        dz, dr = parent_b[1]-parent_a[1], parent_b[2]-parent_a[2]
        length_rz = hypot(dz,dr)
        pieces = max(1,ceil(Int,length_rz/max_boundary_length))
        step_z, step_r = dz/pieces, dr/pieces
        for k in 0:pieces-1
            a = (parent_a[1]+k*step_z,parent_a[2]+k*step_r)
            b = (a[1]+step_z,a[2]+step_r)
            # Lateral area of the conical frustum swept out by this R-Z piece.
            area = pi*(a[2]+b[2])*length_rz/pieces
            area <= 10eps(Float64) && continue
            push!(segments,BoundarySegment(a,b,label,condition,area))
        end
    end
    isempty(segments) && throw(ArgumentError(
        "geometry has no nonzero-area physical boundaries"))
    sort!(segments;by=s -> (s.label,s.a,s.b)) # Deterministic matrix ordering.
    segments
end

function _make_mesh(geometry::AxisymmetricGeometry, boundaries, options::SolverOptions)
    _validate_options(options)
    bc = _validate_boundaries(geometry,boundaries)

    # Give the constrained triangulator a closed boundary-node cycle. Keeping
    # ghost triangles until iteration lets its `each_solid_*` accessors cleanly
    # distinguish the physical domain from the exterior.
    points = copy(geometry.points)
    boundary_nodes = [collect(1:length(points)); 1]
    tri = triangulate(points; boundary_nodes, randomise=false, delete_ghosts=false)
    domain_area = abs(_signed_area(geometry.points))
    max_area = options.max_area > 0 ? options.max_area : domain_area / 200
    refine!(tri; min_angle=options.min_angle, max_area,
            max_points=options.max_mesh_points)

    # DelaunayTriangulation may keep unused points after refinement, so cells—
    # rather than the point array—define which vertices belong to the solution.
    mesh_points = NTuple{2,Float64}[(Float64(p[1]), Float64(p[2])) for p in get_points(tri)]
    cells = NTuple{3,Int}[]
    for T in each_solid_triangle(tri)
        t = (Int(T[1]), Int(T[2]), Int(T[3]))
        a,b,c = mesh_points[t[1]], mesh_points[t[2]], mesh_points[t[3]]
        # Normalize all cells to counter-clockwise orientation. Point location
        # later relies on all three edge-orientation tests being nonnegative.
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
        # Pappus's centroid theorem: V = A * 2πr_centroid. This is exact for a
        # linear triangular R-Z cell, including cells adjacent to the axis.
        push!(volumes, 2pi * center[2] * _triangle_area(a,b,c))
    end

    # Transport segments are generated independently of the triangulator. This
    # prevents volume refinement from multiplying the boundary-exchange matrix.
    # The automatic length tracks the cell scale; specifying a positive value
    # fully decouples boundary and volume resolution.
    boundary_length = options.max_boundary_length > 0 ?
                      options.max_boundary_length : 2sqrt(max_area)
    segments = _make_boundary_segments(geometry,bc,boundary_length)
    return RZMesh(mesh_points, cells, centers, volumes, segments)
end
