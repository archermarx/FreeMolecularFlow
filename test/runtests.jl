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
end

@testset "uniform reservoir equilibrium" begin
    geometry = box_geometry(["axis","reservoir","reservoir","reservoir"])
    pressure, temperature = 0.01, 300.0
    result = solve(geometry,
        Dict("axis"=>Axis(),"reservoir"=>BackPressure(pressure,temperature)),
        XENON;options=QUICK)
    exact_density = pressure/(KB*temperature)
    @test all(isfinite,result.density)
    @test maximum(abs.(result.density ./ exact_density .- 1)) < 0.08
    @test maximum(abs.(result.velocity)) / mean_molecular_speed(XENON,temperature) < 0.03
    @test result.particle_balance_residual < 1e-9
    @test result.exchange_closure_error < 1e-9
    @test maximum(abs.(result.direct_view_factors["reservoir"] .- 1)) < 0.08
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

@testset "configuration and VTK" begin
    mktempdir() do dir
        config = joinpath(@__DIR__,"..","examples","hall_channel.toml")
        loaded = load_config(config)
        @test loaded.gas.molecular_mass ≈ XENON.molecular_mass
        @test loaded.boundaries["anode"] isa Inflow

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
    end
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
