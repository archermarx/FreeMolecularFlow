const Vec3 = NTuple{3,Float64}

# Tiny tuple-based vector operations avoid allocating temporary Vector objects
# in the innermost visibility and solid-angle loops.
_vadd(a::Vec3,b::Vec3) = (a[1]+b[1],a[2]+b[2],a[3]+b[3])
_vsub(a::Vec3,b::Vec3) = (a[1]-b[1],a[2]-b[2],a[3]-b[3])
_vscale(a::Vec3,s::Real) = (a[1]*s,a[2]*s,a[3]*s)
_dot(a::Vec3,b::Vec3) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
_cross(a::Vec3,b::Vec3) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
_norm(a::Vec3) = sqrt(_dot(a,a))
_unit(a::Vec3) = _vscale(a, inv(_norm(a)))

struct SurfacePatch
    vertices::NTuple{3,Vec3}
    center::Vec3
    normal::Vec3
    area::Float64       # Exact-area-rescaled integration weight [m²].
    geometric_area::Float64 # Area of the flat triangulated approximation [m²].
    segment::Int        # Parent BoundarySegment index.
end

# Nodes with left == 0 are leaves and carry patch indices. Interior nodes store
# child indices only. Array indices replace pointers for compact traversal.
struct BVHNode
    lo::Vec3
    hi::Vec3
    left::Int
    right::Int
    patches::Vector{Int}
end

struct PatchBVH
    nodes::Vector{BVHNode}
    root::Int
    patches::Vector{SurfacePatch}
end

# Embed (z,r) at azimuth θ in Cartesian (z,x,y). Keeping z first means output
# velocity components naturally remain ordered (axial, radial, azimuthal).
_xyz(p, theta) = (p[1], p[2]*cos(theta), p[2]*sin(theta))

function _raw_patch(vertices, desired_normal, segment)
    a,b,c = vertices
    cr = _cross(_vsub(b,a), _vsub(c,a))
    area = _norm(cr)/2
    # Wedges meeting the symmetry axis can collapse one of their two triangles.
    # Such zero-area pieces carry no flux and must not enter the BVH.
    area <= 100eps(Float64) && return nothing
    normal = _vscale(cr, inv(2area))
    _dot(normal, desired_normal) < 0 && (normal = _vscale(normal,-1))
    center = _vscale(_vadd(_vadd(a,b),c), 1/3)
    SurfacePatch(vertices, center, normal, area, area, segment)
end

function _revolve_segments(segments::Vector{BoundarySegment}, ntheta::Int)
    # Each R-Z segment sweeps out a conical frustum. Approximate that frustum by
    # two triangles per azimuthal wedge for visibility tests and solid angles.
    patches = SurfacePatch[]
    for (is,s) in pairs(segments)
        first_patch = length(patches)+1
        dz, dr = s.b[1]-s.a[1], s.b[2]-s.a[2]
        len = hypot(dz,dr)
        for k in 0:ntheta-1
            t0, t1 = 2pi*k/ntheta, 2pi*(k+1)/ntheta
            tm = (t0+t1)/2
            # CCW polygon has its domain on the left of each original edge.
            desired = (-dr/len, (dz/len)*cos(tm), (dz/len)*sin(tm))
            a0,b0,b1,a1 = _xyz(s.a,t0),_xyz(s.b,t0),_xyz(s.b,t1),_xyz(s.a,t1)
            for verts in ((a0,b0,b1),(a0,b1,a1))
                patch = _raw_patch(verts, desired, is)
                patch === nothing || push!(patches, patch)
            end
        end
        last_patch = length(patches)
        # A polygonal ring has slightly less area than the smooth surface of
        # revolution. Rescale integration weights to the analytic frustum area;
        # geometry remains flat only for ray intersection purposes.
        approximate_area = sum(patches[i].geometric_area for i in first_patch:last_patch)
        scale = s.area / approximate_area
        for i in first_patch:last_patch
            p = patches[i]
            patches[i] = SurfacePatch(p.vertices,p.center,p.normal,
                                      p.geometric_area*scale,p.geometric_area,p.segment)
        end
    end
    patches
end

function _patch_bounds(p::SurfacePatch)
    lo = ntuple(k -> minimum(v[k] for v in p.vertices), 3)
    hi = ntuple(k -> maximum(v[k] for v in p.vertices), 3)
    lo, hi
end

function _build_bvh(patches::Vector{SurfacePatch}; leaf_size=8)
    # Median-split the longest bounding-box axis. This is inexpensive to build
    # and substantially reduces the O(N) triangle tests for every sight line.
    nodes = BVHNode[]
    function build(indices)
        lows_highs = (_patch_bounds(patches[i]) for i in indices)
        bounds = collect(lows_highs)
        lo = ntuple(k -> minimum(b[1][k] for b in bounds), 3)
        hi = ntuple(k -> maximum(b[2][k] for b in bounds), 3)
        slot = length(nodes)+1
        push!(nodes, BVHNode(lo,hi,0,0,Int[]))
        if length(indices) <= leaf_size
            nodes[slot] = BVHNode(lo,hi,0,0,collect(indices))
        else
            extent = _vsub(hi,lo)
            axis = argmax(extent)
            sorted = sort(collect(indices); by=i -> patches[i].center[axis])
            mid = length(sorted) ÷ 2
            left = build(view(sorted,1:mid))
            right = build(view(sorted,mid+1:length(sorted)))
            nodes[slot] = BVHNode(lo,hi,left,right,Int[])
        end
        slot
    end
    PatchBVH(nodes, build(eachindex(patches)), patches)
end

function _segment_hits_box(origin::Vec3, direction::Vec3, lo::Vec3, hi::Vec3, tol)
    # Slab intersection over t ∈ [0,1]. `direction` is the complete displacement
    # from origin to endpoint rather than a unit ray.
    tmin, tmax = 0.0, 1.0
    for k in 1:3
        if abs(direction[k]) < eps(Float64)
            (origin[k] < lo[k]-tol || origin[k] > hi[k]+tol) && return false
        else
            invd = inv(direction[k])
            t1, t2 = (lo[k]-origin[k])*invd, (hi[k]-origin[k])*invd
            t1 > t2 && ((t1,t2) = (t2,t1))
            tmin, tmax = max(tmin,t1), min(tmax,t2)
            tmin > tmax && return false
        end
    end
    true
end

function _ray_triangle_t(origin::Vec3, direction::Vec3, patch::SurfacePatch, tol)
    # Two-sided Möller–Trumbore intersection. Visibility cares whether any wall
    # blocks the segment, independent of which side of that wall is encountered.
    a,b,c = patch.vertices
    edge1, edge2 = _vsub(b,a), _vsub(c,a)
    h = _cross(direction,edge2)
    det = _dot(edge1,h)
    abs(det) <= tol && return nothing
    f = inv(det)
    s = _vsub(origin,a)
    u = f*_dot(s,h)
    (-tol <= u <= 1+tol) || return nothing
    q = _cross(s,edge1)
    v = f*_dot(direction,q)
    (v >= -tol && u+v <= 1+tol) || return nothing
    f*_dot(edge2,q)
end

function _occluded(bvh::PatchBVH, origin::Vec3, endpoint::Vec3,
                   ignore1::Int, ignore2::Int, tol::Float64)
    direction = _vsub(endpoint,origin)
    stack = Int[bvh.root]
    while !isempty(stack)
        node = bvh.nodes[pop!(stack)]
        _segment_hits_box(origin,direction,node.lo,node.hi,tol) || continue
        if node.left == 0
            for ip in node.patches
                (ip == ignore1 || ip == ignore2) && continue
                t = _ray_triangle_t(origin,direction,bvh.patches[ip],tol)
                # Ignore intersections at either endpoint; these are normally
                # the emitting/receiving patches themselves or shared edges.
                t === nothing || (tol < t < 1-tol && return true)
            end
        else
            push!(stack,node.left,node.right)
        end
    end
    false
end

function _solid_angle(point::Vec3, patch::SurfacePatch)
    # Van Oosterom–Strackee's stable closed form for a triangular solid angle.
    # atan(y,x) retains the correct quadrant for large apparent triangles.
    a,b,c = (_vsub(v,point) for v in patch.vertices)
    la,lb,lc = _norm(a),_norm(b),_norm(c)
    numerator = abs(_dot(a,_cross(b,c)))
    denominator = la*lb*lc + _dot(a,b)*lc + _dot(b,c)*la + _dot(c,a)*lb
    abs(2atan(numerator,denominator))
end

function _vector_solid_angle(point::Vec3, patch::SurfacePatch)
    # Integral of the unit direction vector over a spherical triangle. Each
    # great-circle edge contributes its unit normal times half its arc angle.
    # This moment yields bulk velocity without a center-ray approximation.
    u = ntuple(i -> _unit(_vsub(patch.vertices[i],point)),3)
    value = (0.0,0.0,0.0)
    for (a,b) in ((u[1],u[2]),(u[2],u[3]),(u[3],u[1]))
        cr = _cross(a,b)
        sine = _norm(cr)
        sine <= 10eps(Float64) && continue
        angle = atan(sine,_dot(a,b))
        value = _vadd(value,_vscale(cr,0.5*angle/sine))
    end
    # Orient toward the viewed patch irrespective of its vertex ordering.
    toward_patch = _unit(_vsub(patch.center,point))
    _dot(value,toward_patch) < 0 && (value = _vscale(value,-1))
    value
end

_area_scale(p::SurfacePatch) = p.area / p.geometric_area

function _boundary_exchange(mesh::RZMesh, patches, bvh, tol, ntheta, reporter=nothing)
    # conductance[i,j] has units of area and equals A_i H_ij. Storing this
    # reciprocal form makes the eventual symmetry condition explicit.
    ns = length(mesh.boundary_segments)
    conductance = zeros(ns,ns)
    wedge = 2pi/ntheta
    # Rotational symmetry means a single receiver wedge represents its whole
    # ring. Emitters still span 2π, retaining the full axisymmetric visibility.
    receivers = Int[]
    for (ip,p) in pairs(patches)
        theta = mod(atan(p.center[3],p.center[2]),2pi)
        theta <= wedge + 100eps(Float64) && push!(receivers,ip)
    end
    for (iteration,ip) in enumerate(receivers)
        p = patches[ip]
        for jp in eachindex(patches)
            ip == jp && continue
            q = patches[jp]
            direction = _vsub(q.center,p.center)
            distance = _norm(direction)
            distance <= tol && continue
            ray = _vscale(direction,inv(distance))
            # Both gas-side normals must face the line connecting the patches.
            # q sees the opposite ray direction, hence its leading minus sign.
            cos_p, cos_q = _dot(p.normal,ray), -_dot(q.normal,ray)
            (cos_p > 0 && cos_q > 0) || continue
            _occluded(bvh,p.center,q.center,ip,jp,tol) && continue
            omega_q = _solid_angle(p.center,q) * _area_scale(q)
            # Diffuse (Lambertian) exchange: dH = cos(θ_receiver)dΩ_emitter/π.
            # The receiver patch represents every identical azimuthal wedge.
            g = ntheta * p.area*cos_p*omega_q/pi
            g > 0 || continue
            conductance[p.segment,q.segment] += g
        end
        _status!(reporter,:boundary_exchange,iteration;total=length(receivers))
    end
    # The two independently integrated directions differ slightly at finite
    # quadrature order. Symmetrisation enforces reciprocity exactly.
    conductance .= (conductance .+ transpose(conductance)) ./ 2
    areas = getfield.(mesh.boundary_segments,:area)
    # The centroid/solid-angle rule under-resolves grazing exchange between
    # adjacent rings. Symmetric matrix balancing restores the exact enclosure
    # rule while retaining non-negativity, zeros, and discrete reciprocity.
    d = ones(length(areas))
    balance_iteration = 0
    error = Inf
    for iteration in 1:10_000
        balance_iteration = iteration
        rows = d .* (conductance*d)
        any(rows .<= 0) && throw(ErrorException(
            "boundary quadrature has an isolated surface; increase azimuthal_divisions"))
        error = maximum(abs.(rows ./ areas .- 1))
        _status!(reporter,:exchange_balance,iteration;
                 exchange_closure=error)
        error < 1e-12 && break
        d .*= sqrt.(areas ./ rows)
    end
    conductance .*= d .* transpose(d)
    H = conductance ./ areas # row i divided by receiving area i
    closure = maximum(abs.(sum(conductance;dims=1)[:] ./ areas .- 1))
    _status!(reporter,:exchange_balance,balance_iteration;
             exchange_closure=closure,force=true)
    H, closure
end
