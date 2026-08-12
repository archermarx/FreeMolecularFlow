function _vtk_names(labels)
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

function _required(table, key, context)
    haskey(table,key) || throw(ArgumentError("missing `$key` in $context"))
    table[key]
end

function _parse_boundary(table, label)
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
    (;geometry,boundaries,gas,options,output)
end

"""Run a TOML configuration and write its VTK output."""
function run_config(path::AbstractString)
    config = load_config(path)
    result = solve(config.geometry,config.boundaries,config.gas;options=config.options)
    output = isabspath(config.output) ? config.output : joinpath(dirname(abspath(path)),config.output)
    mkpath(dirname(output))
    written = write_vtk(output,result)
    @printf("cells: %d\n",length(result.mesh.cells))
    @printf("radiosity residual: %.6e\n",result.radiosity_residual)
    @printf("particle balance residual: %.6e\n",result.particle_balance_residual)
    @printf("exchange closure error: %.6e\n",result.exchange_closure_error)
    println("wrote: $written")
    result
end
