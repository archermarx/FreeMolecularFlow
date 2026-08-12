#!/usr/bin/env julia

using FreeMolecularFlow

function usage(io=stdout)
    println(io,"usage: julia --project=. bin/free_molecular_flow.jl [--status-interval SECONDS] CONFIG.toml")
    println(io,"       --status-interval=SECONDS is also accepted; zero disables status output")
end

function parse_args(args)
    interval = 0.0
    config = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h","--help")
            usage()
            exit(0)
        elseif arg == "--status-interval"
            i == length(args) && throw(ArgumentError("--status-interval requires a value"))
            i += 1
            interval = parse(Float64,args[i])
        elseif startswith(arg,"--status-interval=")
            interval = parse(Float64,split(arg,"=";limit=2)[2])
        elseif startswith(arg,"-")
            throw(ArgumentError("unknown option `$arg`"))
        elseif config === nothing
            config = arg
        else
            throw(ArgumentError("only one CONFIG.toml may be supplied"))
        end
        i += 1
    end
    config === nothing && throw(ArgumentError("CONFIG.toml is required"))
    interval >= 0 || throw(ArgumentError("--status-interval must be nonnegative"))
    config, interval
end

try
    config, interval = parse_args(ARGS)
    run_config(config;status_interval=interval)
catch err
    showerror(stderr,err)
    println(stderr)
    usage(stderr)
    exit(1)
end
