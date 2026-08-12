# Package precompilation normally records only methods reached while loading the
# module. A tiny physical solve exercises the deeply nested triangulation and
# visibility specializations that would otherwise dominate the first CLI run.
@setup_workload begin
    geometry = AxisymmetricGeometry(
        [(0.0,0.0),(0.02,0.0),(0.02,0.01),(0.0,0.01)],
        ["axis","outlet","wall","inlet"])
    boundaries = Dict{String,BoundaryCondition}(
        "axis" => Axis(),
        "outlet" => BackPressure(0.0,300.0),
        "wall" => DiffuseWall(400.0),
        "inlet" => Inflow(1e-7,500.0))
    gas = Gas(131.293;unit=:amu)
    options = SolverOptions(
        max_area=2e-4,max_boundary_length=0.02,
        azimuthal_divisions=4,max_mesh_points=100)

    @compile_workload begin
        prepared = prepare(geometry,boundaries;options)
        solve(prepared,gas)
        mktemp() do path,io
            close(io)
            save_prepared(path,prepared)
            load_prepared(path,geometry;options)
        end
    end
end
