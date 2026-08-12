const PREPARED_CACHE_MAGIC = "FreeMolecularFlow.PreparedSolver"
const PREPARED_CACHE_VERSION = 1

_package_version() = Base.pkgversion(@__MODULE__)

function _cache_header()
    (;magic=PREPARED_CACHE_MAGIC,
      format_version=PREPARED_CACHE_VERSION,
      julia_version=VERSION,
      package_version=_package_version())
end

function _validate_cache_header(header)
    header isa NamedTuple &&
        get(header,:magic,nothing) == PREPARED_CACHE_MAGIC ||
        throw(ArgumentError("not a FreeMolecularFlow prepared-solver cache"))
    header.format_version == PREPARED_CACHE_VERSION || throw(ArgumentError(
        "unsupported prepared-solver cache format $(header.format_version); " *
        "expected $(PREPARED_CACHE_VERSION)"))
    header.julia_version == VERSION || throw(ArgumentError(
        "prepared-solver cache was written by Julia $(header.julia_version), " *
        "but this process uses Julia $(VERSION)"))
    header.package_version == _package_version() || throw(ArgumentError(
        "prepared-solver cache was written by FreeMolecularFlow " *
        "$(header.package_version), but this package is $(_package_version())"))
end

function _same_options(a::SolverOptions,b::SolverOptions)
    all(name -> isequal(getfield(a,name),getfield(b,name)),fieldnames(SolverOptions))
end

function _validate_cached_problem(prepared,geometry,options)
    prepared.geometry.points == geometry.points &&
        prepared.geometry.edge_labels == geometry.edge_labels ||
        throw(ArgumentError("prepared-solver cache geometry does not match the requested geometry"))
    _same_options(prepared.options,options) || throw(ArgumentError(
        "prepared-solver cache options do not match the requested SolverOptions"))
    prepared
end


function _cached_problem_matches(prepared,geometry,options)
    prepared.geometry.points == geometry.points &&
        prepared.geometry.edge_labels == geometry.edge_labels &&
        _same_options(prepared.options,options)
end

"""Atomically save reusable geometry-dependent solver data to disk."""
function save_prepared(filename::AbstractString,prepared::PreparedSolver)
    destination = abspath(filename)
    mkpath(dirname(destination))
    temporary,io = mktemp(dirname(destination))
    try
        serialize(io,_cache_header())
        serialize(io,prepared)
        close(io)
        # POSIX rename replaces an existing file atomically, so readers see the
        # complete old cache or the complete new one, never a partial payload.
        Base.Filesystem.rename(temporary,destination)
    catch
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary;force=true)
        rethrow()
    end
    destination
end

"""Load a prepared solver, rejecting incompatible cache or package versions."""
function load_prepared(filename::AbstractString)
    try
        open(filename,"r") do io
            _validate_cache_header(deserialize(io))
            prepared = deserialize(io)
            prepared isa PreparedSolver || throw(ArgumentError(
                "cache payload is not a PreparedSolver"))
            prepared
        end
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError(
            "could not read prepared-solver cache `$(filename)`: $(sprint(showerror,error))"))
    end
end

function load_prepared(filename::AbstractString,geometry::AxisymmetricGeometry;
                       options::SolverOptions=SolverOptions())
    _validate_cached_problem(load_prepared(filename),geometry,options)
end

"""
    prepare_cached(filename, geometry, boundaries; options, ...)

Load a matching `PreparedSolver` from `filename`, or prepare and atomically
write one when the file is absent or stale. Corrupt caches are reported rather
than silently discarded.
"""
function prepare_cached(filename::AbstractString,
                        geometry::AxisymmetricGeometry,boundaries;
                        options::SolverOptions=SolverOptions(),kwargs...)
    if isfile(filename)
        prepared = load_prepared(filename)
        _cached_problem_matches(prepared,geometry,options) && return prepared
    end
    prepared = prepare(geometry,boundaries;options,kwargs...)
    save_prepared(filename,prepared)
    prepared
end
