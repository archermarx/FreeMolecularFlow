#!/usr/bin/env julia

using FreeMolecularFlow

if length(ARGS) != 1
    println(stderr,"usage: julia --project=. bin/free_molecular_flow.jl CONFIG.toml")
    exit(2)
end

try
    run_config(ARGS[1])
catch err
    showerror(stderr,err,catch_backtrace())
    println(stderr)
    exit(1)
end
