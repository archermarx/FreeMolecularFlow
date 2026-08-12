function _vtk_names(labels)
    # Boundary labels are user text and may contain spaces or punctuation.
    # Produce stable, collision-free identifiers suitable for VTK and CSV.
    used = Set{String}()
    names = Dict{String,String}()
    for label in labels
        base = replace(lowercase(strip(label)), r"[^a-z0-9_]+" => "_")
        isempty(base) && (base = "boundary")
        name, suffix = base, 2
        while name in used
            name = "$(base)_$(suffix)"
            suffix += 1
        end
        push!(used,name)
        names[label] = name
    end
    names
end

"""Write a `FlowResult` as an R-Z unstructured VTK XML (`.vtu`) file."""
function write_vtk(filename::AbstractString, result::FlowResult)
    # Write the meridional mesh at y=0. ParaView still expects 3-D coordinates
    # and vectors, so points and velocities retain three components.
    root = endswith(lowercase(filename),".vtu") ? filename[1:end-4] : String(filename)
    points = zeros(3,length(result.mesh.points))
    for (i,p) in pairs(result.mesh.points)
        points[1,i], points[2,i] = p
    end
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE,collect(t)) for t in result.mesh.cells]
    names = _vtk_names(result.labels)
    files = vtk_grid(root,points,cells) do vtk
        vtk["number_density",VTKCellData()] = result.density
        vtk["velocity",VTKCellData()] = result.velocity
        for label in result.labels
            safe = names[label]
            vtk["direct_view_factor_$(safe)",VTKCellData()] = result.direct_view_factors[label]
            vtk["density_from_$(safe)",VTKCellData()] = result.density_contributions[label]
        end
        vtk["radiosity_residual",VTKFieldData()] = [result.radiosity_residual]
        vtk["particle_balance_residual",VTKFieldData()] = [result.particle_balance_residual]
        vtk["exchange_closure_error",VTKFieldData()] = [result.exchange_closure_error]
        vtk["azimuthal_divisions",VTKFieldData()] = [result.azimuthal_divisions]
        vtk["azimuthal_convergence_error",VTKFieldData()] =
            [result.azimuthal_convergence_error]
    end
    only(files)
end

struct CellLocator
    mesh::RZMesh
    lo::NTuple{2,Float64}
    hi::NTuple{2,Float64}
    dimensions::NTuple{2,Int}
    bin_scale::NTuple{2,Float64}
    bins::Vector{Vector{Int}}
    tolerance::Float64
end

_bin_index(value,origin,scale,count) =
    clamp(floor(Int,(value-origin)*scale)+1,1,count)

function _cell_contains(mesh::RZMesh,cell::Int,point,tolerance)
    t = mesh.cells[cell]
    a,b,c = mesh.points[t[1]],mesh.points[t[2]],mesh.points[t[3]]
    o1,o2,o3 = _orient(a,b,point),_orient(b,c,point),_orient(c,a,point)
    o1 >= -tolerance && o2 >= -tolerance && o3 >= -tolerance
end

function _build_cell_locator(mesh::RZMesh)
    # A compact uniform grid gives nearly constant-time point candidates for
    # large extraction paths. Triangle bounding boxes are inserted into every
    # bin they overlap, so points on cell or bin edges remain robust.
    scale = maximum(abs,Iterators.flatten(mesh.points);init=1.0)
    tolerance = 1e-12*max(scale,1.0)^2
    lo = (minimum(first,mesh.points),minimum(last,mesh.points))
    hi = (maximum(first,mesh.points),maximum(last,mesh.points))
    width = max(hi[1]-lo[1],eps(Float64))
    height = max(hi[2]-lo[2],eps(Float64))
    aspect = width/height
    nx = max(1,round(Int,sqrt(length(mesh.cells)*aspect)))
    ny = max(1,ceil(Int,length(mesh.cells)/nx))
    bin_scale = (nx/width,ny/height)
    bins = [Int[] for _ in 1:nx*ny]
    for (cell,t) in pairs(mesh.cells)
        a,b,c = mesh.points[t[1]],mesh.points[t[2]],mesh.points[t[3]]
        ix0 = _bin_index(min(a[1],b[1],c[1]),lo[1],bin_scale[1],nx)
        ix1 = _bin_index(max(a[1],b[1],c[1]),lo[1],bin_scale[1],nx)
        iy0 = _bin_index(min(a[2],b[2],c[2]),lo[2],bin_scale[2],ny)
        iy1 = _bin_index(max(a[2],b[2],c[2]),lo[2],bin_scale[2],ny)
        for iy in iy0:iy1, ix in ix0:ix1
            push!(bins[ix+nx*(iy-1)],cell)
        end
    end
    CellLocator(mesh,lo,hi,(nx,ny),bin_scale,bins,tolerance)
end

function _containing_cell(locator::CellLocator,point::NTuple{2,Float64})
    lo,hi = locator.lo,locator.hi
    tol = locator.tolerance
    (lo[1]-tol <= point[1] <= hi[1]+tol &&
     lo[2]-tol <= point[2] <= hi[2]+tol) || return 0
    nx,ny = locator.dimensions
    ix = _bin_index(point[1],lo[1],locator.bin_scale[1],nx)
    iy = _bin_index(point[2],lo[2],locator.bin_scale[2],ny)
    for cell in locator.bins[ix+nx*(iy-1)]
        _cell_contains(locator.mesh,cell,point,tol) && return cell
    end
    0
end

function _csv_value(io::IO, value::Real)
    # Sixteen significant digits preserve Float64 data through a CSV round trip.
    # Non-finite values deliberately become empty fields rather than nonstandard
    # NaN/Inf spellings that downstream spreadsheet tools interpret differently.
    isfinite(value) ? @printf(io,"%.16g",value) : nothing
end

function _path_samples(line::ExtractionLine)
    # Parameterize the entire polyline by cumulative arc length. This makes both
    # `num_points` and `spacing` uniform across corners rather than restarting
    # the requested resolution independently on every segment.
    segment_lengths = [hypot(line.points[i+1][1]-line.points[i][1],
                             line.points[i+1][2]-line.points[i][2])
                       for i in 1:length(line.points)-1]
    cumulative = [0.0; cumsum(segment_lengths)]
    total = last(cumulative)
    distances = if line.num_points !== nothing
        collect(range(0.0,total;length=line.num_points))
    else
        values = collect(0.0:line.spacing:total)
        # Always include the final control point. Unless spacing divides the
        # length exactly, the final interval is intentionally shorter.
        if total-last(values) > 10eps(total)
            push!(values,total)
        else
            values[end] = total
        end
        values
    end
    samples = NamedTuple[]
    for (sample,distance) in pairs(distances)
        # At an exact corner, associate the sample with the segment ending at
        # that corner. Clamp the final endpoint to the final valid segment.
        segment = min(searchsortedlast(cumulative,distance),length(segment_lengths))
        local_fraction = (distance-cumulative[segment])/segment_lengths[segment]
        a,b = line.points[segment],line.points[segment+1]
        point = (a[1]+local_fraction*(b[1]-a[1]),
                 a[2]+local_fraction*(b[2]-a[2]))
        push!(samples,(;sample,segment,fraction=distance/total,distance,point))
    end
    samples
end

function _extraction_columns(line::ExtractionLine, labels, names)
    columns = ["sample","path_segment","fraction","distance","z","r",
               "inside_domain","cell_index"]
    :number_density in line.fields && push!(columns,"number_density")
    :velocity in line.fields && append!(columns,["velocity_z","velocity_r","velocity_theta"])
    if :view_factors in line.fields
        append!(columns,["direct_view_factor_$(names[label])" for label in labels])
    end
    if :density_contributions in line.fields
        append!(columns,["density_from_$(names[label])" for label in labels])
    end
    columns
end

function _write_extracted_values(io,line,labels,values,index)
    if :number_density in line.fields
        print(io,','); _csv_value(io,values.density[index])
    end
    if :velocity in line.fields
        for component in 1:3
            print(io,','); _csv_value(io,values.velocity[component,index])
        end
    end
    if :view_factors in line.fields
        for label in labels
            print(io,','); _csv_value(io,values.direct[label][index])
        end
    end
    if :density_contributions in line.fields
        for label in labels
            print(io,','); _csv_value(io,values.contributions[label][index])
        end
    end
end

"""
    write_extraction_line(filename, result, line)

Uniformly sample a straight or piecewise-linear path and write selected fields
to CSV. Direct evaluation uses the Katz view-factor moments at the exact sample
locations; cell evaluation uses the containing triangle's cell-centered value.
"""
function write_extraction_line(filename::AbstractString, result::FlowResult,
                               line::ExtractionLine; evaluator=nothing,
                               locator=nothing,reporter=nothing)
    names = _vtk_names(result.labels)
    columns = _extraction_columns(line,result.labels,names)
    samples = _path_samples(line)
    # Domain membership is determined from the R-Z solution mesh even in direct
    # mode. Direct evaluation changes field accuracy, not what counts as gas.
    locator === nothing && (locator = _build_cell_locator(result.mesh))
    cells = [_containing_cell(locator,sample.point) for sample in samples]
    outside = findall(==(0),cells)
    if line.outside_domain === :error && !isempty(outside)
        first_outside = samples[first(outside)]
        throw(ArgumentError("extraction line `$(line.name)` leaves the gas domain " *
            "at sample $(first_outside.sample), point $(first_outside.point)"))
    end
    kept = line.outside_domain === :drop ? findall(!=(0),cells) : collect(eachindex(samples))

    inside_indices = findall(!=(0),cells)
    cell_values = (;density=result.density,velocity=result.velocity,
                   direct=result.direct_view_factors,
                   contributions=result.density_contributions)
    direct_values = nothing
    direct_lookup = zeros(Int,length(samples))
    _status!(reporter,:line_extraction,0;total=length(inside_indices),force=true)
    if line.method === :direct && !isempty(inside_indices)
        # Evaluate only interior points; outside samples, when retained, receive
        # empty solution columns and never participate in expensive ray tracing.
        evaluator === nothing && (evaluator = result.evaluator)
        points = [samples[i].point for i in inside_indices]
        _,direct,contributions,density,velocity = _evaluate_fields(
            result.mesh,points,evaluator.patches,evaluator.occluder,
            result.boundary_flux,result.gas,evaluator.tol;
            ntheta=evaluator.azimuthal_divisions,
            reporter,phase=:line_extraction,
            closure=result.exchange_closure_error,
            radiosity=result.radiosity_residual)
        direct_values = (;direct,contributions,density,velocity)
        # Direct arrays are compacted to interior points, so map original path
        # sample indices back to their row in those arrays during CSV emission.
        direct_lookup[inside_indices] .= eachindex(inside_indices)
    end

    mkpath(dirname(abspath(filename)))
    open(filename,"w") do io
        println(io,join(columns,','))
        for sample_index in kept
            sample = samples[sample_index]
            cell = cells[sample_index]
            print(io,sample.sample,',',sample.segment,','); _csv_value(io,sample.fraction)
            print(io,','); _csv_value(io,sample.distance)
            print(io,','); _csv_value(io,sample.point[1])
            print(io,','); _csv_value(io,sample.point[2])
            print(io,',',cell > 0,',',cell)
            if cell > 0
                values = line.method === :cell ? cell_values : direct_values
                value_index = line.method === :cell ? cell : direct_lookup[sample_index]
                _write_extracted_values(io,line,result.labels,values,value_index)
            else
                for _ in 1:(length(columns)-8)
                    print(io,',')
                end
            end
            println(io)
        end
    end
    line.method === :cell && _status!(reporter,:line_extraction,
        length(inside_indices);total=length(inside_indices),
        exchange_closure=result.exchange_closure_error,
        radiosity=result.radiosity_residual,force=true)
    String(filename)
end

function _parse_extraction_point(table, context)
    table isa AbstractDict || throw(ArgumentError("$context must be {z=..., r=...}"))
    (Float64(_required(table,"z",context)),Float64(_required(table,"r",context)))
end

function _required(table, key, context)
    haskey(table,key) || throw(ArgumentError("missing `$key` in $context"))
    table[key]
end

function _parse_boundary(table, label)
    # A few unambiguous aliases make hand-written TOML forgiving. Specular is
    # recognized only to provide a precise unsupported-physics diagnostic.
    context = "boundaries.$label"
    required(key) = _required(table,key,context)
    kind = lowercase(String(required("type")))
    if kind == "inflow"
        Inflow(required("mass_flow_rate"),required("temperature"))
    elseif kind in ("back_pressure","backpressure","outlet")
        BackPressure(required("pressure"),required("temperature"))
    elseif kind in ("diffuse_wall","diffuse","wall")
        DiffuseWall(required("temperature"))
    elseif kind == "axis"
        Axis()
    elseif kind in ("specular","specular_wall")
        throw(ArgumentError("specular reflection is not implemented; use a diffuse_wall or revise the model"))
    else
        throw(ArgumentError("unknown boundary type `$kind` for label `$label`"))
    end
end

function _parse_gas(table)
    haskey(table,"molecular_mass_amu") && return Gas(table["molecular_mass_amu"];unit=:amu)
    haskey(table,"molecular_mass_kg") && return Gas(table["molecular_mass_kg"])
    throw(ArgumentError("gas requires molecular_mass_amu or molecular_mass_kg"))
end

_case_path(config_path,path) =
    isabspath(path) ? path : joinpath(dirname(abspath(config_path)),path)

"""Load geometry, boundary conditions, gas, solver options, and output path from TOML."""
function load_config(path::AbstractString)
    # Parse into validated domain types immediately. The solver therefore never
    # needs to branch on raw TOML dictionaries or perform unit conversion.
    cfg = TOML.parsefile(path)
    gtab = _required(cfg,"geometry","configuration")
    points = _required(gtab,"points","geometry")
    labels = _required(gtab,"edge_labels","geometry")
    geometry = AxisymmetricGeometry(points,labels)

    gas = _parse_gas(_required(cfg,"gas","configuration"))

    btab = _required(cfg,"boundaries","configuration")
    boundaries = Dict{String,BoundaryCondition}(
        String(label) => _parse_boundary(table,label) for (label,table) in btab)

    opt = get(cfg,"solver",Dict{String,Any}())
    options = SolverOptions(
        max_area=Float64(get(opt,"max_area",0.0)),
        max_boundary_length=Float64(get(opt,"max_boundary_length",0.0)),
        min_angle=Float64(get(opt,"min_angle",20.0)),
        azimuthal_divisions=Int(get(opt,"azimuthal_divisions",64)),
        azimuthal_tolerance=Float64(get(opt,"azimuthal_tolerance",0.0)),
        max_azimuthal_divisions=Int(get(opt,"max_azimuthal_divisions",256)),
        radiosity_tolerance=Float64(get(opt,"radiosity_tolerance",1e-10)),
        visibility_tolerance=Float64(get(opt,"visibility_tolerance",1e-10)),
        max_mesh_points=Int(get(opt,"max_mesh_points",200_000)))
    output = String(get(get(cfg,"output",Dict{String,Any}()),"path","free_molecular_flow.vtu"))
    cache_table = get(cfg,"cache",nothing)
    cache = cache_table === nothing ? nothing :
            String(_required(cache_table,"path","cache"))
    # TOML [[extraction_lines]] entries arrive as a vector of dictionaries.
    # Named coordinate tables avoid the easy-to-miss positional [z,r] ordering.
    extraction_lines = ExtractionLine[]
    for (i,table) in enumerate(get(cfg,"extraction_lines",Any[]))
        context = "extraction_lines[$i]"
        point_tables = _required(table,"points",context)
        points = [_parse_extraction_point(p,"$context.points[$j]")
                  for (j,p) in pairs(point_tables)]
        push!(extraction_lines,ExtractionLine(_required(table,"name",context),points;
            num_points=get(table,"num_points",nothing),
            spacing=get(table,"spacing",nothing),
            method=get(table,"method","direct"),
            outside_domain=get(table,"outside_domain","keep"),
            fields=get(table,"fields",EXTRACTION_FIELDS),
            filename=get(table,"filename",nothing)))
    end
    names = getfield.(extraction_lines,:name)
    length(unique(names)) == length(names) ||
        throw(ArgumentError("extraction-line names must be unique"))
    filenames = getfield.(extraction_lines,:filename)
    length(unique(filenames)) == length(filenames) ||
        throw(ArgumentError("extraction-line filenames must be unique"))
    (;geometry,boundaries,gas,options,output,cache,extraction_lines)
end

"""Run a TOML configuration and write its VTK output."""
function run_config(path::AbstractString; status_interval::Real=0.0,
                    status_io::IO=stdout)
    config = load_config(path)
    reporter = _make_reporter(status_interval,status_io,nothing)
    cache_path = config.cache === nothing ? nothing : _case_path(path,config.cache)
    result = if cache_path === nothing
        solve(config.geometry,config.boundaries,config.gas;
              options=config.options,status_reporter=reporter)
    else
        prepared = prepare_cached(cache_path,config.geometry,config.boundaries;
                                  options=config.options,status_reporter=reporter)
        solve(prepared,config.boundaries,config.gas;status_reporter=reporter)
    end
    # All relative outputs are case-relative, not process-working-directory
    # relative. A case can therefore be launched reliably from any directory.
    output = _case_path(path,config.output)
    mkpath(dirname(output))
    written = write_vtk(output,result)
    extraction_files = String[]
    # Surface construction is independent of the extraction path; share it
    # whenever more than one line requests exact direct evaluation.
    evaluator = any(line -> line.method === :direct,config.extraction_lines) ?
                result.evaluator : nothing
    locator = isempty(config.extraction_lines) ? nothing :
              _build_cell_locator(result.mesh)
    for line in config.extraction_lines
        path_out = _case_path(path,line.filename)
        push!(extraction_files,write_extraction_line(path_out,result,line;
                                                     evaluator,locator,reporter))
    end
    @printf("cells: %d\n",length(result.mesh.cells))
    @printf("radiosity residual: %.6e\n",result.radiosity_residual)
    @printf("particle balance residual: %.6e\n",result.particle_balance_residual)
    @printf("exchange closure error: %.6e\n",result.exchange_closure_error)
    @printf("azimuthal divisions: %d\n",result.azimuthal_divisions)
    isfinite(result.azimuthal_convergence_error) && @printf(
        "azimuthal convergence error: %.6e\n",result.azimuthal_convergence_error)
    cache_path === nothing || println("cache: $cache_path")
    println("wrote: $written")
    for file in extraction_files
        println("wrote: $file")
    end
    result
end
