const Vec3 = NTuple{3,Float64}

# Tiny tuple-based vector operations avoid allocating temporary Vector objects
# in the innermost visibility and solid-angle loops.
_vadd(a::Vec3,b::Vec3) = (a[1]+b[1],a[2]+b[2],a[3]+b[3])
_vsub(a::Vec3,b::Vec3) = (a[1]-b[1],a[2]-b[2],a[3]-b[3])
_vscale(a::Vec3,s::Real) = (a[1]*s,a[2]*s,a[3]*s)
_dot(a::Vec3,b::Vec3) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
_cross(a::Vec3,b::Vec3) = (a[2]*b[3]-a[3]*b[2],
                           a[3]*b[1]-a[1]*b[3],
                           a[1]*b[2]-a[2]*b[1])
_norm(a::Vec3) = sqrt(_dot(a,a))

struct SurfacePatch
    origin::Vec3         # First triangle vertex.
    edge1::Vec3          # Triangle edge b-a.
    edge2::Vec3          # Triangle edge c-a.
    center::Vec3
    surface_center::Vec3 # Center projected onto the smooth revolved segment.
    normal::Vec3
    area::Float64       # Exact-area-rescaled integration weight [m²].
    area_scale::Float64 # Ratio of smooth frustum area to polygonal area.
    segment::Int        # Parent BoundarySegment index.
end

_patch_vertices(p::SurfacePatch) =
    (p.origin,_vadd(p.origin,p.edge1),_vadd(p.origin,p.edge2))

# Nodes with left == 0 are leaves and reference a contiguous range in the BVH's
# shared patch-index array. Avoiding a separately allocated Vector per leaf
# improves locality during the millions of visibility traversals.
struct BVHNode
    lo::Vec3
    hi::Vec3
    left::Int
    right::Int
    first_patch::Int
    patch_count::Int
end

struct PatchBVH
    nodes::Vector{BVHNode}
    root::Int
    patches::Vector{SurfacePatch}
    patch_indices::Vector{Int}
end

const CONICAL_SURFACE = UInt8(1)
const DISK_SURFACE = UInt8(2)

"""One original R-Z polygon edge revolved about the z axis."""
struct RevolvedSurface
    kind::UInt8
    z0::Float64
    r0::Float64
    slope::Float64
    zmin::Float64
    zmax::Float64
    rmin::Float64
    rmax::Float64
    lo::Vec3
    hi::Vec3
end

struct AnalyticOccluder
    nodes::Vector{BVHNode}
    root::Int
    surfaces::Vector{RevolvedSurface}
    surface_indices::Vector{Int}
end

# Values shared by every bounding-box and triangle test along one sight line.
# In particular, direction reciprocals belong to the ray—not to each BVH node.
struct RaySegment
    origin::Vec3
    direction::Vec3
    inverse_direction::Vec3
    parallel_axes::NTuple{3,Bool}
end

@inline function _ray_segment(origin::Vec3,endpoint::Vec3)
    direction = _vsub(endpoint,origin)
    parallel = ntuple(i -> abs(direction[i]) < eps(Float64),3)
    inverse = ntuple(i -> parallel[i] ? 0.0 : inv(direction[i]),3)
    RaySegment(origin,direction,inverse,parallel)
end

# Embed (z,r) at azimuth θ in Cartesian (z,x,y). Keeping z first means output
_xyz(p, theta) = (p[1], p[2]*cos(theta), p[2]*sin(theta))

function _smooth_surface_center(center,s::BoundarySegment)
    dz = s.b[1]-s.a[1]
    abs(dz) <= eps(Float64) && return center # Annular disk.
    radius = s.a[2] + (center[1]-s.a[1])*(s.b[2]-s.a[2])/dz
    polygonal_radius = hypot(center[2],center[3])
    polygonal_radius <= eps(Float64) && return center
    scale = radius/polygonal_radius
    (center[1],center[2]*scale,center[3]*scale)
end


function _smooth_surface_normal(center,s::BoundarySegment)
    dz,dr = s.b[1]-s.a[1],s.b[2]-s.a[2]
    length_rz = hypot(dz,dr)
    radius = hypot(center[2],center[3])
    radial_x = radius <= eps(Float64) ? 1.0 : center[2]/radius
    radial_y = radius <= eps(Float64) ? 0.0 : center[3]/radius
    (-dr/length_rz,(dz/length_rz)*radial_x,(dz/length_rz)*radial_y)
end

function _raw_patch(vertices,desired_normal,segment,s::BoundarySegment)
    a,b,c = vertices
    edge1,edge2 = _vsub(b,a),_vsub(c,a)
    cr = _cross(edge1,edge2)
    area = _norm(cr)/2
    # Wedges meeting the symmetry axis can collapse one of their two triangles.
    # Such zero-area pieces carry no flux and must not enter the BVH.
    area <= 100eps(Float64) && return nothing
    center = _vscale(_vadd(_vadd(a,b),c), 1/3)
    smooth_center = _smooth_surface_center(center,s)
    normal = _smooth_surface_normal(smooth_center,s)
    # Desired and analytic normals should agree; retain this orientation guard
    # for nearly degenerate apex triangles.
    _dot(normal,desired_normal) < 0 && (normal = _vscale(normal,-1))
    SurfacePatch(a,edge1,edge2,center,smooth_center,normal,area,1.0,segment)
end

function _revolve_segments(segments::Vector{BoundarySegment}, ntheta::Int)
    # Each R-Z segment sweeps out a conical frustum. Triangles provide solid-angle
    # integration elements; visibility uses the exact smooth surface separately.
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
            # For an even number of wedges, mirror the diagonal across the
            # meridional y=0 plane in the negative-y half. This makes every
            # positive-y triangle have an exact reflected partner, allowing
            # field reconstruction to integrate only half of each ring.
            mirrored_half = iseven(ntheta) && k >= ntheta ÷ 2
            triangles = mirrored_half ?
                ((a0,b0,a1),(a1,b0,b1)) :
                ((a0,b0,b1),(a0,b1,a1))
            for verts in triangles
                patch = _raw_patch(verts,desired,is,s)
                patch === nothing || push!(patches, patch)
            end
        end
        last_patch = length(patches)
        # A polygonal ring has slightly less area than the smooth surface of
        # revolution. Rescale integration weights to the analytic frustum area;
        # geometry remains flat only for ray intersection purposes.
        approximate_area = sum(patches[i].area for i in first_patch:last_patch)
        scale = s.area / approximate_area
        for i in first_patch:last_patch
            p = patches[i]
            patches[i] = SurfacePatch(
                p.origin,p.edge1,p.edge2,p.center,p.surface_center,p.normal,
                p.area*scale,scale,p.segment)
        end
    end
    patches
end

function _patch_bounds(p::SurfacePatch)
    vertices = _patch_vertices(p)
    lo = ntuple(k -> minimum(v[k] for v in vertices),3)
    hi = ntuple(k -> maximum(v[k] for v in vertices),3)
    lo, hi
end

function _build_bvh(patches::Vector{SurfacePatch}; leaf_size=8)
    # Retained as a triangle-visibility reference for diagnostics and regression
    # tests. Production solves use the much smaller analytic surface hierarchy.
    nodes = BVHNode[]
    patch_indices = Int[]
    function build(indices)
        lows_highs = (_patch_bounds(patches[i]) for i in indices)
        bounds = collect(lows_highs)
        lo = ntuple(k -> minimum(b[1][k] for b in bounds), 3)
        hi = ntuple(k -> maximum(b[2][k] for b in bounds), 3)
        slot = length(nodes)+1
        push!(nodes,BVHNode(lo,hi,0,0,0,0))
        if length(indices) <= leaf_size
            first_patch = length(patch_indices)+1
            append!(patch_indices,indices)
            nodes[slot] = BVHNode(lo,hi,0,0,first_patch,length(indices))
        else
            extent = _vsub(hi,lo)
            axis = argmax(extent)
            sorted = sort(collect(indices); by=i -> patches[i].center[axis])
            mid = length(sorted) ÷ 2
            left = build(view(sorted,1:mid))
            right = build(view(sorted,mid+1:length(sorted)))
            nodes[slot] = BVHNode(lo,hi,left,right,0,0)
        end
        slot
    end
    root = build(eachindex(patches))
    PatchBVH(nodes,root,patches,patch_indices)
end

function _revolved_surfaces(geometry::AxisymmetricGeometry)
    surfaces = RevolvedSurface[]
    for i in eachindex(geometry.points)
        a = geometry.points[i]
        b = geometry.points[mod1(i+1,length(geometry.points))]
        a[2] == 0 && b[2] == 0 && continue # The symmetry axis has zero area.
        dz = b[1]-a[1]
        zmin,zmax = minmax(a[1],b[1])
        rmin,rmax = minmax(a[2],b[2])
        kind = abs(dz) <= eps(Float64) ? DISK_SURFACE : CONICAL_SURFACE
        slope = kind == DISK_SURFACE ? 0.0 : (b[2]-a[2])/dz
        lo = (zmin,-rmax,-rmax)
        hi = (zmax,rmax,rmax)
        push!(surfaces,RevolvedSurface(
            kind,a[1],a[2],slope,zmin,zmax,rmin,rmax,lo,hi))
    end
    isempty(surfaces) && throw(ArgumentError(
        "geometry has no nonzero-area surfaces for visibility testing"))
    surfaces
end

function _build_occluder(geometry::AxisymmetricGeometry; leaf_size=2)
    surfaces = _revolved_surfaces(geometry)
    nodes = BVHNode[]
    surface_indices = Int[]
    function build(indices)
        lo = ntuple(k -> minimum(surfaces[i].lo[k] for i in indices),3)
        hi = ntuple(k -> maximum(surfaces[i].hi[k] for i in indices),3)
        slot = length(nodes)+1
        push!(nodes,BVHNode(lo,hi,0,0,0,0))
        if length(indices) <= leaf_size
            first_surface = length(surface_indices)+1
            append!(surface_indices,indices)
            nodes[slot] = BVHNode(
                lo,hi,0,0,first_surface,length(indices))
        else
            # Full surfaces of revolution all have x/y centroid zero, so z is
            # the only useful partition axis for this very small hierarchy.
            sorted = sort(collect(indices);
                by=i -> (surfaces[i].zmin+surfaces[i].zmax)/2)
            mid = length(sorted) ÷ 2
            left = build(view(sorted,1:mid))
            right = build(view(sorted,mid+1:length(sorted)))
            nodes[slot] = BVHNode(lo,hi,left,right,0,0)
        end
        slot
    end
    root = build(eachindex(surfaces))
    AnalyticOccluder(nodes,root,surfaces,surface_indices)
end

@inline function _segment_box_entry(ray::RaySegment,lo::Vec3,hi::Vec3,tol)
    # Slab intersection over t ∈ [0,1]. `direction` is the complete displacement
    # from origin to endpoint rather than a unit ray.
    tmin, tmax = 0.0, 1.0
    for k in 1:3
        if ray.parallel_axes[k]
            (ray.origin[k] < lo[k]-tol || ray.origin[k] > hi[k]+tol) && return Inf
        else
            t1 = (lo[k]-ray.origin[k])*ray.inverse_direction[k]
            t2 = (hi[k]-ray.origin[k])*ray.inverse_direction[k]
            t1 > t2 && ((t1,t2) = (t2,t1))
            tmin, tmax = max(tmin,t1), min(tmax,t2)
            tmin > tmax && return Inf
        end
    end
    tmin
end

@inline function _ray_triangle_t(ray::RaySegment,patch::SurfacePatch,tol)
    # Two-sided Möller–Trumbore intersection. Visibility cares whether any wall
    # blocks the segment, independent of which side of that wall is encountered.
    a = patch.origin
    edge1,edge2 = patch.edge1,patch.edge2
    h = _cross(ray.direction,edge2)
    det = _dot(edge1,h)
    abs(det) <= tol && return nothing
    f = inv(det)
    s = _vsub(ray.origin,a)
    u = f*_dot(s,h)
    (-tol <= u <= 1+tol) || return nothing
    q = _cross(s,edge1)
    v = f*_dot(ray.direction,q)
    (v >= -tol && u+v <= 1+tol) || return nothing
    f*_dot(edge2,q)
end

function _new_traversal_stack(bvh)
    # The median-split tree needs only logarithmic depth; 64 entries already
    # cover physically unrealizable patch counts. `_push_stack!` still grows the
    # array defensively if the BVH construction strategy ever changes.
    Vector{Int}(undef,min(max(length(bvh.nodes),1),64))
end

@inline function _push_stack!(stack,top,node)
    top += 1
    top > length(stack) && resize!(stack,2length(stack))
    @inbounds stack[top] = node
    top
end

function _occluded(bvh::PatchBVH, origin::Vec3, endpoint::Vec3,
                   ignore1::Int, ignore2::Int, tol::Float64,
                   stack::Vector{Int})
    ray = _ray_segment(origin,endpoint)
    root = bvh.nodes[bvh.root]
    isfinite(_segment_box_entry(ray,root.lo,root.hi,tol)) ||
        return false
    top = _push_stack!(stack,0,bvh.root)
    while top > 0
        @inbounds node_index = stack[top]
        top -= 1
        node = bvh.nodes[node_index]
        if node.left == 0
            last_patch = node.first_patch+node.patch_count-1
            for slot in node.first_patch:last_patch
                ip = bvh.patch_indices[slot]
                (ip == ignore1 || ip == ignore2) && continue
                t = _ray_triangle_t(ray,bvh.patches[ip],tol)
                # Ignore intersections at either endpoint; these are normally
                # the emitting/receiving patches themselves or shared edges.
                t === nothing || (tol < t < 1-tol && return true)
            end
        else
            left,right = bvh.nodes[node.left],bvh.nodes[node.right]
            left_entry = _segment_box_entry(ray,left.lo,left.hi,tol)
            right_entry = _segment_box_entry(ray,right.lo,right.hi,tol)
            # LIFO traversal pushes the farther child first. When an occluder is
            # present this tends to find the closest blocking patch sooner.
            if left_entry <= right_entry
                isfinite(right_entry) && (top = _push_stack!(stack,top,node.right))
                isfinite(left_entry) && (top = _push_stack!(stack,top,node.left))
            else
                isfinite(left_entry) && (top = _push_stack!(stack,top,node.left))
                isfinite(right_entry) && (top = _push_stack!(stack,top,node.right))
            end
        end
    end
    false
end

@inline _interior_intersection(t,tol) = isfinite(t) && tol < t < 1-tol

@inline function _surface_z_contains(surface,ray,t,tol)
    z = ray.origin[1] + t*ray.direction[1]
    surface.zmin-tol <= z <= surface.zmax+tol
end

@inline function _hits_disk(surface,ray,tol)
    abs(ray.direction[1]) <= eps(Float64) && return false
    t = (surface.z0-ray.origin[1])/ray.direction[1]
    _interior_intersection(t,tol) || return false
    x = ray.origin[2] + t*ray.direction[2]
    y = ray.origin[3] + t*ray.direction[3]
    radius = hypot(x,y)
    surface.rmin-tol <= radius <= surface.rmax+tol
end

@inline function _cone_root_hits(surface,ray,t,tol)
    _interior_intersection(t,tol) &&
        _surface_z_contains(surface,ray,t,tol)
end

@inline function _hits_cone(surface,ray,tol)
    # rho(t)^2 = r(z(t))^2 gives a quadratic for finite cones and cylinders.
    z,x,y = ray.origin
    dz,dx,dy = ray.direction
    radius0 = surface.r0 + surface.slope*(z-surface.z0)
    radius_step = surface.slope*dz
    A = dx*dx + dy*dy - radius_step*radius_step
    B = 2(x*dx + y*dy - radius0*radius_step)
    C = x*x + y*y - radius0*radius0
    scale = max(abs(A),abs(B),abs(C),floatmin(Float64))
    threshold = 64eps(Float64)*scale
    if abs(A) <= threshold
        abs(B) <= threshold && return false # Parallel to or contained in surface.
        return _cone_root_hits(surface,ray,-C/B,tol)
    end

    discriminant = B*B - 4A*C
    discriminant_scale = max(abs(B*B),abs(4A*C),floatmin(Float64))
    discriminant < -64eps(Float64)*discriminant_scale && return false
    root = sqrt(max(discriminant,0.0))
    # This form avoids cancellation when B and sqrt(discriminant) are close.
    q = -0.5*(B+copysign(root,B))
    t1 = q/A
    _cone_root_hits(surface,ray,t1,tol) && return true
    t2 = q == 0 ? -B/(2A) : C/q
    _cone_root_hits(surface,ray,t2,tol)
end

@inline function _hits_surface(surface,ray,tol)
    surface.kind == DISK_SURFACE ? _hits_disk(surface,ray,tol) :
                                  _hits_cone(surface,ray,tol)
end

function _occluded(bvh::AnalyticOccluder,origin::Vec3,endpoint::Vec3,
                   ignore1::Int,ignore2::Int,tol::Float64,
                   stack::Vector{Int})
    ray = _ray_segment(origin,endpoint)
    root = bvh.nodes[bvh.root]
    isfinite(_segment_box_entry(ray,root.lo,root.hi,tol)) || return false
    top = _push_stack!(stack,0,bvh.root)
    while top > 0
        @inbounds node_index = stack[top]
        top -= 1
        node = bvh.nodes[node_index]
        if node.left == 0
            last_surface = node.first_patch+node.patch_count-1
            for slot in node.first_patch:last_surface
                surface = bvh.surfaces[bvh.surface_indices[slot]]
                _hits_surface(surface,ray,tol) && return true
            end
        else
            left,right = bvh.nodes[node.left],bvh.nodes[node.right]
            left_entry = _segment_box_entry(ray,left.lo,left.hi,tol)
            right_entry = _segment_box_entry(ray,right.lo,right.hi,tol)
            if left_entry <= right_entry
                isfinite(right_entry) && (top = _push_stack!(stack,top,node.right))
                isfinite(left_entry) && (top = _push_stack!(stack,top,node.left))
            else
                isfinite(left_entry) && (top = _push_stack!(stack,top,node.left))
                isfinite(right_entry) && (top = _push_stack!(stack,top,node.right))
            end
        end
    end
    false
end

# Retain a convenient allocation-owning method for infrequent diagnostic use.
# Performance-critical loops pass reusable scratch storage explicitly.
function _occluded(bvh::PatchBVH, origin::Vec3, endpoint::Vec3,
                   ignore1::Int, ignore2::Int, tol::Float64)
    _occluded(bvh,origin,endpoint,ignore1,ignore2,tol,
              _new_traversal_stack(bvh))
end

function _occluded(bvh::AnalyticOccluder,origin::Vec3,endpoint::Vec3,
                   ignore1::Int,ignore2::Int,tol::Float64)
    _occluded(bvh,origin,endpoint,ignore1,ignore2,tol,
              _new_traversal_stack(bvh))
end

@inline function _patch_rays(point::Vec3, patch::SurfacePatch)
    rays = map(vertex -> _vsub(vertex,point),_patch_vertices(patch))
    rays, map(_norm,rays)
end

@inline function _scalar_solid_angle(rays, lengths)
    a,b,c = rays
    la,lb,lc = lengths
    numerator = abs(_dot(a,_cross(b,c)))
    denominator = la*lb*lc + _dot(a,b)*lc + _dot(b,c)*la + _dot(c,a)*lb
    # Two-argument atan retains the correct quadrant for large apparent triangles.
    abs(2atan(numerator,denominator))
end

function _solid_angle(point::Vec3, patch::SurfacePatch)
    _scalar_solid_angle(_patch_rays(point,patch)...)
end

function _solid_angle_moments(point::Vec3, patch::SurfacePatch)
    # This is deliberately fused: field reconstruction needs the scalar, vector,
    # and diagonal second moments, and sharing vertex rays and lengths is
    # measurably faster in this hot loop.
    a,b,c = (_vsub(vertex,point) for vertex in _patch_vertices(patch))
    la,lb,lc = _norm(a),_norm(b),_norm(c)
    numerator = abs(_dot(a,_cross(b,c)))
    denominator = la*lb*lc + _dot(a,b)*lc + _dot(b,c)*la + _dot(c,a)*lb
    omega = abs(2atan(numerator,denominator))

    directions = (_vscale(a,inv(la)),_vscale(b,inv(lb)),_vscale(c,inv(lc)))
    vector_omega = (0.0,0.0,0.0)
    # On the unit sphere, Δ_s(s_i*s_j) = 2δ_ij - 6s_i*s_j. Applying the
    # surface-divergence theorem converts the second solid-angle moment into
    # exact great-circle edge integrals. Only the diagonal is retained because
    # these are the three directional temperatures requested by the R-Z model.
    second_omega = (omega/3,omega/3,omega/3)
    for (u,v,w) in ((directions[1],directions[2],directions[3]),
                    (directions[2],directions[3],directions[1]),
                    (directions[3],directions[1],directions[2]))
        cross_uv = _cross(u,v)
        sine = _norm(cross_uv)
        sine <= 10eps(Float64) && continue
        angle = atan(sine,_dot(u,v))
        vector_omega = _vadd(
            vector_omega,_vscale(cross_uv,0.5*angle/sine))

        # Choose the edge conormal that points into the spherical triangle.
        inward = _vscale(cross_uv,inv(sine))
        _dot(inward,w) < 0 && (inward = _vscale(inward,-1))
        # Integral of the unit direction along the minor great-circle arc. This
        # midpoint form remains stable when the edge angle approaches π.
        midpoint = _vadd(u,v)
        midpoint_norm = _norm(midpoint)
        midpoint_norm <= 10eps(Float64) && continue
        edge_integral = _vscale(
            midpoint,2sin(angle/2)/midpoint_norm)
        second_omega = ntuple(i ->
            second_omega[i] + inward[i]*edge_integral[i]/3,3)
    end
    _dot(vector_omega,_vsub(patch.center,point)) < 0 &&
        (vector_omega = _vscale(vector_omega,-1))
    omega,vector_omega,second_omega
end

function _receiver_patches(patches,ns,ntheta)
    wedge = 2pi/ntheta
    receivers = [Int[] for _ in 1:ns]
    for (ip,p) in pairs(patches)
        theta = mod(atan(p.center[3],p.center[2]),2pi)
        theta <= wedge+100eps(Float64) && push!(receivers[p.segment],ip)
    end
    receivers
end

function _patch_conductance(receiver,emitter,bvh,receiver_index,emitter_index,
                            ntheta,tol,tolerance_squared,stack)
    direction = _vsub(emitter.surface_center,receiver.surface_center)
    distance_squared = _dot(direction,direction)
    distance_squared > tolerance_squared || return 0.0

    # Test the two half-spaces before paying for a square root or ray traversal.
    receiver_projection = _dot(receiver.normal,direction)
    emitter_projection = -_dot(emitter.normal,direction)
    (receiver_projection > 0 && emitter_projection > 0) || return 0.0
    _occluded(bvh,receiver.surface_center,emitter.surface_center,
              receiver_index,emitter_index,
              tol,stack) && return 0.0

    cosine = receiver_projection/sqrt(distance_squared)
    omega = _solid_angle(receiver.surface_center,emitter)*emitter.area_scale
    # The receiver patch represents all identical azimuthal wedges.
    max(ntheta*receiver.area*cosine*omega/pi,0.0)
end

function _balance_conductance!(conductance,areas,reporter)
    # Finite quadrature makes the two directions differ slightly. Symmetrize
    # for reciprocity, then scale rows/columns together to close the enclosure.
    conductance .= (conductance .+ transpose(conductance)) ./ 2
    scaling = ones(length(areas))
    iteration = 0
    for current in 1:10_000
        iteration = current
        row_sums = scaling .* (conductance*scaling)
        any(row_sums .<= 0) && throw(ErrorException(
            "boundary quadrature has an isolated surface; increase azimuthal_divisions"))
        error = maximum(abs.(row_sums ./ areas .- 1))
        _status!(reporter,:exchange_balance,current;exchange_closure=error)
        error < 1e-12 && break
        scaling .*= sqrt.(areas ./ row_sums)
    end
    conductance .*= scaling .* transpose(scaling)
    closure = maximum(abs.(vec(sum(conductance;dims=1)) ./ areas .- 1))
    _status!(reporter,:exchange_balance,iteration;
             exchange_closure=closure,force=true)
    conductance ./ areas, closure
end

function _boundary_exchange(mesh::RZMesh,patches,bvh,tol,ntheta,reporter=nothing)
    # conductance[i,j] has units of area and equals A_i H_ij. Storing this
    # reciprocal form makes the eventual symmetry condition explicit.
    ns = length(mesh.boundary_segments)
    conductance = zeros(ns,ns)
    # Rotational symmetry means a single receiver wedge represents its whole
    # ring. Emitters still span 2π, retaining the full axisymmetric visibility.
    receivers = _receiver_patches(patches,ns,ntheta)
    # A segment owns one output row, so receiver groups can run concurrently
    # without atomics or thread-local conductance matrices. Accumulation order
    # within each row remains deterministic across thread counts.
    traversal_stacks = [_new_traversal_stack(bvh) for _ in 1:Threads.maxthreadid()]
    completed = Threads.Atomic{Int}(0)
    reporter_lock = ReentrantLock()
    tolerance_squared = tol^2
    Threads.@threads for iseg in 1:ns
        traversal_stack = traversal_stacks[Threads.threadid()]
        for ip in receivers[iseg]
            receiver = patches[ip]
            for jp in eachindex(patches)
                ip == jp && continue
                emitter = patches[jp]
                g = _patch_conductance(receiver,emitter,bvh,ip,jp,ntheta,tol,
                                       tolerance_squared,traversal_stack)
                conductance[iseg,emitter.segment] += g
            end
        end
        _thread_status!(reporter,reporter_lock,completed,:boundary_exchange,ns)
    end
    areas = getfield.(mesh.boundary_segments,:area)
    _balance_conductance!(conductance,areas,reporter)
end
