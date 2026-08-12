module FreeMolecularFlow

# External packages are kept deliberately small: triangulation creates the
# meridional mesh, while WriteVTK handles the only binary output format.
using DelaunayTriangulation
using LinearAlgebra
using PrecompileTools
using Printf
using Serialization
using TOML
using WriteVTK

export AxisymmetricGeometry, Gas, Inflow, BackPressure, DiffuseWall, Axis,
       SolverOptions, PreparedSolver, FlowResult, prepare, solve,
       prepare_cached, save_prepared, load_prepared,
       write_vtk, load_config, run_config,
       ExtractionLine, write_extraction_line, number_density,
       mean_molecular_speed

# Include order follows the dependency chain. In particular, visibility.jl
# operates on the mesh types constructed by geometry.jl, and solver.jl combines
# both layers before io.jl exposes configuration and output functions.
include("types.jl")       # Public data model and small geometric predicates.
include("geometry.jl")    # R-Z polygon validation and triangulation.
include("visibility.jl")  # Revolved surfaces, ray tracing, and view factors.
include("solver.jl")      # Radiosity and field-moment reconstruction.
include("cache.jl")       # Versioned cross-process PreparedSolver reuse.
include("io.jl")          # TOML, VTK, CSV, and command-line orchestration.
include("precompile.jl")  # Representative workload for low-latency CLI starts.

end
