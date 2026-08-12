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
    end
    only(files)
end

function _containing_cell(mesh::RZMesh, point::NTuple{2,Float64})
    # Cells were normalized counter-clockwise during meshing. A point lies in a
    # triangle when it is on the left of (or numerically on) every directed edge.
    # Returning zero gives callers a cheap sentinel for outside-domain samples.
    scale = maximum(abs,Iterators.flatten(mesh.points);init=1.0)
    tolerance = 1e-12*max(scale,1.0)^2
    for (i,t) in pairs(mesh.cells)
        a,b,c = mesh.points[t[1]],mesh.points[t[2]],mesh.points[t[3]]
        o1,o2,o3 = _orient(a,b,point),_orient(b,c,point),_orient(c,a,point)
        o1 >= -tolerance && o2 >= -tolerance && o3 >= -tolerance && return i
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

function _prepare_direct_evaluator(result::FlowResult)
    # FlowResult stores the refined boundary and solver options, but not the
    # relatively bulky patch/BVH cache. Rebuild it once and share it across all
    # direct extraction paths requested by a configuration.
    patches = _revolve_segments(result.mesh.boundary_segments,
                                result.options.azimuthal_divisions)
    bvh = _build_bvh(patches)
    scale = maximum(abs,Iterators.flatten(result.mesh.points);init=1.0)
    tol = result.options.visibility_tolerance*max(scale,1.0)
    (;patches,bvh,tol)
end

function _write_extracted_values(io, line, labels, density, velocity,
                                 direct, contributions, index)
    if :number_density in line.fields
        print(io,','); _csv_value(io,density[index])
    end
    if :velocity in line.fields
        for component in 1:3
            print(io,','); _csv_value(io,velocity[component,index])
        end
    end
    if :view_factors in line.fields
        for label in labels
            print(io,','); _csv_value(io,direct[label][index])
        end
    end
    if :density_contributions in line.fields
        for label in labels
            print(io,','); _csv_value(io,contributions[label][index])
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
                               reporter=nothing)
    names = _vtk_names(result.labels)
    columns = _extraction_columns(line,result.labels,names)
    samples = _path_samples(line)
    # Domain membership is determined from the R-Z solution mesh even in direct
    # mode. Direct evaluation changes field accuracy, not what counts as gas.
    cells = [_containing_cell(result.mesh,sample.point) for sample in samples]
    outside = findall(==(0),cells)
    if line.outside_domain === :error && !isempty(outside)
        first_outside = samples[first(outside)]
        throw(ArgumentError("extraction line `$(line.name)` leaves the gas domain " *
            "at sample $(first_outside.sample), point $(first_outside.point)"))
    end
    kept = line.outside_domain === :drop ? findall(!=(0),cells) : collect(eachindex(samples))

    inside_indices = findall(!=(0),cells)
    direct_values = nothing
    direct_lookup = Dict{Int,Int}()
    _status!(reporter,:line_extraction,0;total=length(inside_indices),force=true)
    if line.method === :direct && !isempty(inside_indices)
        # Evaluate only interior points; outside samples, when retained, receive
        # empty solution columns and never participate in expensive ray tracing.
        evaluator === nothing && (evaluator = _prepare_direct_evaluator(result))
        points = [samples[i].point for i in inside_indices]
        _,direct,contributions,density,velocity = _evaluate_fields(
            result.mesh,points,evaluator.patches,evaluator.bvh,
            result.boundary_flux,result.gas,evaluator.tol;
            reporter,phase=:line_extraction,
            closure=result.exchange_closure_error,
            radiosity=result.radiosity_residual)
        direct_values = (;direct,contributions,density,velocity)
        # Direct arrays are compacted to interior points, so map original path
        # sample indices back to their row in those arrays during CSV emission.
        direct_lookup = Dict(sample_index=>value_index
                             for (value_index,sample_index) in pairs(inside_indices))
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
                if line.method === :cell
                    _write_extracted_values(io,line,result.labels,result.density,
                        result.velocity,result.direct_view_factors,
                        result.density_contributions,cell)
                else
                    value_index = direct_lookup[sample_index]
                    _write_extracted_values(io,line,result.labels,
                        direct_values.density,direct_values.velocity,
                        direct_values.direct,direct_values.contributions,value_index)
                end
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
    kind = lowercase(String(_required(table,"type","boundaries.$label")))
    if kind == "inflow"
        Inflow(_required(table,"mass_flow_rate","boundaries.$label"),
               _required(table,"temperature","boundaries.$label"))
    elseif kind in ("back_pressure","backpressure","outlet")
        BackPressure(_required(table,"pressure","boundaries.$label"),
                     _required(table,"temperature","boundaries.$label"))
    elseif kind in ("diffuse_wall","diffuse","wall")
        DiffuseWall(_required(table,"temperature","boundaries.$label"))
    elseif kind == "axis"
        Axis()
    elseif kind in ("specular","specular_wall")
        throw(ArgumentError("specular reflection is not implemented; use a diffuse_wall or revise the model"))
    else
        throw(ArgumentError("unknown boundary type `$kind` for label `$label`"))
    end
end

"""Load geometry, boundary conditions, gas, solver options, and output path from TOML."""
function load_config(path::AbstractString)
    # Parse into validated domain types immediately. The solver therefore never
    # needs to branch on raw TOML dictionaries or perform unit conversion.
    cfg = TOML.parsefile(path)
    gtab = _required(cfg,"geometry","configuration")
    points = _required(gtab,"points","geometry")
    labels = _required(gtab,"edge_labels","geometry")
    geometry = AxisymmetricGeometry(points,labels)

    gastab = _required(cfg,"gas","configuration")
    if haskey(gastab,"molecular_mass_amu")
        gas = Gas(gastab["molecular_mass_amu"];unit=:amu)
    elseif haskey(gastab,"molecular_mass_kg")
        gas = Gas(gastab["molecular_mass_kg"])
    else
        throw(ArgumentError("gas requires molecular_mass_amu or molecular_mass_kg"))
    end

    btab = _required(cfg,"boundaries","configuration")
    boundaries = Dict{String,BoundaryCondition}(
        String(label) => _parse_boundary(table,label) for (label,table) in btab)

    opt = get(cfg,"solver",Dict{String,Any}())
    options = SolverOptions(
        max_area=Float64(get(opt,"max_area",0.0)),
        min_angle=Float64(get(opt,"min_angle",20.0)),
        azimuthal_divisions=Int(get(opt,"azimuthal_divisions",64)),
        radiosity_tolerance=Float64(get(opt,"radiosity_tolerance",1e-10)),
        visibility_tolerance=Float64(get(opt,"visibility_tolerance",1e-10)),
        max_mesh_points=Int(get(opt,"max_mesh_points",200_000)))
    output = String(get(get(cfg,"output",Dict{String,Any}()),"path","free_molecular_flow.vtu"))
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
    (;geometry,boundaries,gas,options,output,extraction_lines)
end

"""Run a TOML configuration and write its VTK output."""
function run_config(path::AbstractString; status_interval::Real=0.0,
                    status_io::IO=stdout)
    status_interval >= 0 || throw(ArgumentError("status_interval must be nonnegative"))
    config = load_config(path)
    reporter = status_interval > 0 ? StatusReporter(status_interval,status_io) : nothing
    result = solve(config.geometry,config.boundaries,config.gas;
                   options=config.options,status_reporter=reporter)
    # All relative outputs are case-relative, not process-working-directory
    # relative. A case can therefore be launched reliably from any directory.
    output = isabspath(config.output) ? config.output : joinpath(dirname(abspath(path)),config.output)
    mkpath(dirname(output))
    written = write_vtk(output,result)
    extraction_files = String[]
    # Surface construction is independent of the extraction path; share it
    # whenever more than one line requests exact direct evaluation.
    evaluator = any(line -> line.method === :direct,config.extraction_lines) ?
                _prepare_direct_evaluator(result) : nothing
    for line in config.extraction_lines
        path_out = isabspath(line.filename) ? line.filename :
                   joinpath(dirname(abspath(path)),line.filename)
        push!(extraction_files,write_extraction_line(path_out,result,line;
                                                     evaluator,reporter))
    end
    @printf("cells: %d\n",length(result.mesh.cells))
    @printf("radiosity residual: %.6e\n",result.radiosity_residual)
    @printf("particle balance residual: %.6e\n",result.particle_balance_residual)
    @printf("exchange closure error: %.6e\n",result.exchange_closure_error)
    println("wrote: $written")
    for file in extraction_files
        println("wrote: $file")
    end
    result
end
