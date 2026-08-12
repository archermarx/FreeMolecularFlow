using FreeMolecularFlow
using Test

const KB = 1.380649e-23

function box_geometry(labels=["axis","right","top","left"])
    AxisymmetricGeometry([(0.0,0.0),(0.1,0.0),(0.1,0.05),(0.0,0.05)],labels)
end

const QUICK = SolverOptions(max_area=1e-3,min_angle=20.0,azimuthal_divisions=8)
const XENON = Gas(131.293;unit=:amu)

@testset "types and geometry validation" begin
    @test XENON.molecular_mass > 0
    @test mean_molecular_speed(XENON,300.0) > 0
    @test_throws ArgumentError Gas(0.0)
    @test_throws ArgumentError Inflow(-1.0,300.0)
    @test_throws ArgumentError DiffuseWall(0.0)
    @test_throws ArgumentError AxisymmetricGeometry(
        [(0.0,0.0),(1.0,1.0),(0.0,1.0),(1.0,0.0)],fill("x",4))
    @test_throws ArgumentError AxisymmetricGeometry(
        [(0.0,-0.1),(1.0,0.0),(0.0,1.0)],fill("x",3))
    @test_throws ArgumentError FreeMolecularFlow._validate_options(
        SolverOptions(max_boundary_length=-1.0))
end

@testset "independent boundary resolution" begin
    geometry = box_geometry(["axis","reservoir","reservoir","reservoir"])
    boundaries = Dict("axis"=>Axis(),
                      "reservoir"=>BackPressure(0.01,300.0))
    coarse = FreeMolecularFlow._make_mesh(geometry,boundaries,
        SolverOptions(max_area=1e-3,max_boundary_length=0.0125,
                      azimuthal_divisions=8))
    fine = FreeMolecularFlow._make_mesh(geometry,boundaries,
        SolverOptions(max_area=1e-4,max_boundary_length=0.0125,
                      azimuthal_divisions=8))
    @test length(fine.cells) > length(coarse.cells)
    @test getfield.(fine.boundary_segments,:a) ==
          getfield.(coarse.boundary_segments,:a)
    @test getfield.(fine.boundary_segments,:b) ==
          getfield.(coarse.boundary_segments,:b)
    @test all(s -> hypot(s.b[1]-s.a[1],s.b[2]-s.a[2]) <= 0.0125+eps(),
              fine.boundary_segments)
end

@testset "uniform reservoir equilibrium" begin
    geometry = box_geometry(["axis","reservoir","reservoir","reservoir"])
    pressure, temperature = 0.01, 300.0
    status = IOBuffer()
    result = solve(geometry,
        Dict("axis"=>Axis(),"reservoir"=>BackPressure(pressure,temperature)),
        XENON;options=QUICK,status_interval=1e9,status_io=status)
    exact_density = pressure/(KB*temperature)
    @test all(isfinite,result.density)
    @test maximum(abs.(result.density ./ exact_density .- 1)) < 0.08
    @test maximum(abs.(result.velocity)) / mean_molecular_speed(XENON,temperature) < 0.03
    @test result.particle_balance_residual < 1e-9
    @test result.exchange_closure_error < 1e-9
    @test maximum(abs.(result.direct_view_factors["reservoir"] .- 1)) < 0.08
    status_text = String(take!(status))
    @test occursin("phase",status_text)
    @test occursin("iteration",status_text)
    @test occursin("elapsed(s)",status_text)
    @test occursin("mesh",status_text)
    @test occursin("complete",status_text)
    @test_throws ArgumentError solve(geometry,
        Dict("axis"=>Axis(),"reservoir"=>BackPressure(pressure,temperature)),
        XENON;options=QUICK,status_interval=-1)
end

@testset "inflow, diffuse wall, and linearity" begin
    geometry = box_geometry()
    function run(mdot)
        solve(geometry,Dict(
            "axis"=>Axis(),
            "right"=>BackPressure(0.0,300.0),
            "top"=>DiffuseWall(400.0),
            "left"=>Inflow(mdot,500.0)),XENON;options=QUICK)
    end
    a, b = run(1e-6), run(2e-6)
    @test a.radiosity_residual < 1e-10
    @test a.particle_balance_residual < 1e-9
    @test all(a.density .> 0)
    @test b.density ≈ 2 .* a.density rtol=2e-12
    summed = reduce(+,values(a.density_contributions))
    @test summed ≈ a.density rtol=5e-15
end

@testset "prepared parameter sweeps" begin
    geometry = box_geometry(["axis","reservoir","reservoir","reservoir"])
    boundaries = Dict("axis"=>Axis(),
                      "reservoir"=>BackPressure(0.01,300.0))
    prepared = prepare(geometry,boundaries;options=QUICK)
    @test size(prepared.direction_moments,1) == 2
    base = solve(prepared,XENON)
    doubled_boundaries = Dict("axis"=>Axis(),
                              "reservoir"=>BackPressure(0.02,300.0))
    doubled = solve(prepared,doubled_boundaries,XENON)
    @test doubled.evaluator.bvh === prepared.bvh
    @test doubled.density ≈ 2 .* base.density rtol=5e-14
    @test doubled.velocity ≈ base.velocity rtol=5e-14 atol=1e-12
    @test_throws ArgumentError solve(prepared,
        Dict("axis"=>Axis()),XENON)
end

@testset "optimized field quadrature" begin
    geometry = box_geometry(["axis","reservoir","reservoir","reservoir"])
    boundaries = Dict("axis"=>Axis(),
                      "reservoir"=>BackPressure(0.01,300.0))
    mesh = FreeMolecularFlow._make_mesh(geometry,boundaries,QUICK)
    patches = FreeMolecularFlow._revolve_segments(
        mesh.boundary_segments,QUICK.azimuthal_divisions)
    bvh = FreeMolecularFlow._build_bvh(patches)
    @test sort(bvh.patch_indices) == collect(eachindex(patches))
    @test all(node -> node.left != 0 || node.patch_count > 0,bvh.nodes)
    q = FreeMolecularFlow._prescribed_fluxes(mesh,XENON)
    tol = QUICK.visibility_tolerance

    indices,weight = FreeMolecularFlow._field_patch_quadrature(
        patches,QUICK.azimuthal_divisions)
    @test length(indices)*2 == length(patches)
    @test weight == 2.0

    # The fused routine must reproduce the scalar formula, and its directional
    # moment cannot be larger than the scalar solid angle.
    point = (mesh.centers[1][1],mesh.centers[1][2],0.0)
    patch = patches[first(indices)]
    omega,vector_omega = FreeMolecularFlow._solid_angle_moments(point,patch)
    @test omega ≈ FreeMolecularFlow._solid_angle(point,patch) rtol=2e-15
    @test FreeMolecularFlow._norm(vector_omega) <= omega

    # Half-ring integration should agree with explicitly summing both members
    # of every mirrored pair. Passing an odd marker selects the full-ring path.
    samples = mesh.centers[1:min(3,end)]
    half = FreeMolecularFlow._evaluate_fields(
        mesh,samples,patches,bvh,q,XENON,tol;
        ntheta=QUICK.azimuthal_divisions)
    full = FreeMolecularFlow._evaluate_fields(
        mesh,samples,patches,bvh,q,XENON,tol;ntheta=1)
    @test half[1] == full[1]
    @test half[4] ≈ full[4] rtol=2e-14
    @test half[5] ≈ full[5] rtol=2e-14 atol=1e-12
    for label in half[1]
        @test half[2][label] ≈ full[2][label] rtol=2e-14
        @test half[3][label] ≈ full[3][label] rtol=2e-14
    end

    odd_patches = FreeMolecularFlow._revolve_segments(mesh.boundary_segments,7)
    odd_indices,odd_weight = FreeMolecularFlow._field_patch_quadrature(odd_patches,7)
    @test collect(odd_indices) == collect(eachindex(odd_patches))
    @test odd_weight == 1.0
end

@testset "configuration and VTK" begin
    mktempdir() do dir
        config = joinpath(@__DIR__,"..","examples","hall_channel.toml")
        loaded = load_config(config)
        @test loaded.gas.molecular_mass ≈ XENON.molecular_mass
        @test loaded.boundaries["anode"] isa Inflow
        @test loaded.options.max_boundary_length == 5e-3

        geometry = box_geometry(["axis","reservoir","reservoir","reservoir"])
        result = solve(geometry,
            Dict("axis"=>Axis(),"reservoir"=>BackPressure(0.01,300.0)),
            XENON;options=QUICK)
        path = write_vtk(joinpath(dir,"result.vtu"),result)
        @test isfile(path)
        xml = read(path,String)
        @test occursin("number_density",xml)
        @test occursin("velocity",xml)
        @test occursin("direct_view_factor_reservoir",xml)
        @test occursin("density_from_reservoir",xml)

        line = ExtractionLine("crossing",
            [(-0.01,0.025),(0.11,0.025)];num_points=5,method=:cell,
            outside_domain=:keep,fields=[:number_density,:velocity],
            filename="line.csv")
        csv_path = write_extraction_line(joinpath(dir,line.filename),result,line)
        rows = readlines(csv_path)
        @test length(rows) == 6
        @test startswith(rows[1],"sample,path_segment,fraction,distance,z,r,inside_domain,cell_index,number_density")
        @test split(rows[2],',')[7:8] == ["false","0"]
        @test split(rows[4],',')[7] == "true"
        @test split(rows[6],',')[7:8] == ["false","0"]

        dropped = ExtractionLine("drop outside",
            [(-0.01,0.025),(0.11,0.025)];num_points=5,method=:cell,
            outside_domain=:drop,fields=[:number_density])
        dropped_rows = readlines(write_extraction_line(
            joinpath(dir,dropped.filename),result,dropped))
        @test length(dropped_rows) == 4 # header plus three interior samples

        strict = ExtractionLine("strict",
            [(-0.01,0.025),(0.05,0.025)];num_points=3,method=:cell,
            outside_domain=:error,fields=[:number_density])
        @test_throws ArgumentError write_extraction_line(
            joinpath(dir,strict.filename),result,strict)

        direct_line = ExtractionLine("direct profile",
            [(0.01,0.025),(0.09,0.025)];spacing=0.02,
            fields=[:number_density],outside_domain=:error)
        direct_path = write_extraction_line(joinpath(dir,direct_line.filename),
                                            result,direct_line)
        direct_rows = readlines(direct_path)
        @test length(direct_rows) == 6
        @test direct_rows[1] == "sample,path_segment,fraction,distance,z,r,inside_domain,cell_index,number_density"
        @test all(row -> split(row,',')[7] == "true",direct_rows[2:end])

        locator = FreeMolecularFlow._build_cell_locator(result.mesh)
        @test all(FreeMolecularFlow._containing_cell(locator,center) == i
                  for (i,center) in pairs(result.mesh.centers))
        @test FreeMolecularFlow._containing_cell(locator,(-1.0,0.0)) == 0
        brute_force(point) = something(findfirst(
            cell -> FreeMolecularFlow._cell_contains(
                result.mesh,cell,point,locator.tolerance),eachindex(result.mesh.cells)),0)
        for z in range(-0.01,0.11;length=13), r in range(0.0,0.06;length=9)
            point = (z,r)
            indexed = FreeMolecularFlow._containing_cell(locator,point)
            @test (indexed == 0) == (brute_force(point) == 0)
        end

        polyline = ExtractionLine("bent path",
            [(0.01,0.01),(0.09,0.01),(0.09,0.04)];spacing=0.04,
            method=:cell,fields=[:number_density])
        polyline_rows = readlines(write_extraction_line(
            joinpath(dir,polyline.filename),result,polyline))
        @test length(polyline_rows) == 5
        @test split(polyline_rows[4],',')[2] == "2"
    end
end

@testset "line configuration" begin
    config = load_config(joinpath(@__DIR__,"..","examples","spt100.toml"))
    @test length(config.extraction_lines) == 1
    line = only(config.extraction_lines)
    # The example's physical path may evolve; verify parser/type invariants
    # without coupling this unit test to particular SPT-100 sample coordinates.
    @test !isempty(line.name)
    @test length(line.points) >= 2
    @test all(p -> all(isfinite,p) && p[2] >= 0,line.points)
    @test (line.num_points === nothing) != (line.spacing === nothing)
    @test line.method in (:direct,:cell)
    @test line.outside_domain in (:keep,:drop,:error)
    @test endswith(line.filename,".csv")
    @test_throws ArgumentError ExtractionLine("bad",[(0,0),(1,1)];num_points=1)
    @test_throws ArgumentError ExtractionLine("bad",[(0,0),(1,1)];
                                              num_points=2,spacing=0.1)
    @test_throws ArgumentError ExtractionLine("bad",[(0,0),(1,1)];
                                              spacing=0.1,method=:unknown)
    @test_throws ArgumentError ExtractionLine("bad",[(0,-1),(1,1)];spacing=0.1)
end

@testset "boundary errors" begin
    geometry = box_geometry()
    @test_throws ArgumentError solve(geometry,Dict(
        "axis"=>Axis(),"right"=>BackPressure(0,300),
        "top"=>DiffuseWall(300)),XENON;options=QUICK)
    closed = AxisymmetricGeometry(
        [(0.0,0.01),(0.1,0.01),(0.1,0.05),(0.0,0.05)],fill("wall",4))
    @test_throws ArgumentError solve(closed,
        Dict("wall"=>DiffuseWall(300)),XENON;options=QUICK)
end
