#!/usr/bin/env julia
# Opt-in acquisition/noise evidence, no fit or production selector.
# julia --startup-file=no --threads=1 --project=. test/diagnose_acquisition_noise.jl --help
include(joinpath(@__DIR__, "lib", "acquisition_noise_diagnostics.jl"))
using .AcquisitionNoiseDiagnostics

function acquisition_help(io=stdout)
    println(io,"""
    Usage (Julia 1.13):
      julia --startup-file=no --threads=1 --project=. test/diagnose_acquisition_noise.jl \
        --file scan.sxm --geometry saved_lobe_features.tsv \
        --config frozen_physical.toml --selected-summary original_summary.tsv \
        --settings config/label_free_exploration.toml --outdir NEW_SEPARATE_DIRECTORY

    All six paths are required. No defaults, fitting, N selection or class labels.
    Geometry must contain exactly one actual saved candidate at N_selected for this
    file: file,N,lobe,amplitude,x_nm,y_nm,sigma_parallel_nm,sigma_perp_nm,axis_x,
    axis_y,origin_x_nm,origin_y_nm,baseline,tilt_x,tilt_y (optional skew_ratio).
    Original selected-summary refined_policy is required for adaptive support.
    [acquisition_noise] in --settings specifies all diagnostic search/noise settings.
    The output directory must not exist. It cannot be reused or overwritten.

    Shift: forward[y,x] is compared with backward[y+dy,x+dx], in analysis-grid
    pixels after frozen preprocessing.stride. Native SXM backward x reversal is
    already applied. Only signed positive correlation is a plausible alignment;
    boundary and ambiguous optima remain unresolved. No automatic range expansion.
    Background estimates exclude native ROI and guarded saved molecular footprints
    at all searched shifts. Fitted-residual diagnostics are not noise estimates.
    ACF/block counts are descriptive, not independent observations or n_eff.
    Explicit auxiliary_channel adds direct channel concordance on the same Z
    support/lag. Z-auxiliary correlations describe coupling, not chemical evidence.
    No [1,Z,dZ/dx,dZ/dy] nuisance projection is performed in this first pass.
    Multifile/heavy real diagnostics belong on Viper, not this local workstation.
    """)
end

function acquisition_options(args)
    required=("--file","--geometry","--config","--selected-summary","--settings","--outdir")
    opts=Dict{String,String}()
    i=1
    while i<=length(args)
        arg=args[i]
        arg in required || error("Unknown acquisition diagnostic argument: $arg")
        haskey(opts,arg) && error("Duplicate argument: $arg")
        i<length(args) || error("Missing value for $arg")
        startswith(args[i+1],"--") && error("Missing value for $arg")
        opts[arg]=args[i+1]
        i+=2
    end
    all(k->haskey(opts,k),required) || error("All six input/output arguments are required; use --help")
    for k in required[1:5]
        isfile(opts[k]) || error("Input is not a file: $(opts[k])")
    end
    out=opts["--outdir"]
    (ispath(out) || islink(out)) && error("Exclusive output directory already exists: $out")
    return opts
end

function acquisition_main(args=ARGS)
    if args in (["--help"],["-h"])
        acquisition_help()
        return nothing
    end
    v"1.13.0" <= VERSION < v"1.14.0" || error("This diagnostic requires native Julia 1.13")
    opts=acquisition_options(args)
    s=load_settings(opts["--settings"])
    case=load_case(opts["--file"],opts["--geometry"],opts["--config"],opts["--selected-summary"],s)
    d=diagnose(case,s)
    inputs=Dict(k[3:end]=>v for (k,v) in opts if k!="--outdir")
    write_diagnostic(opts["--outdir"],case,d,s;inputs=inputs)
    println("Diagnostic only: ",opts["--outdir"])
    println("Registration: ",d.registration.status,"; constant pixels=",d.registration.n)
    println("Background: ",d.noise_fwd.status," / ",d.noise_bwd.status)
    println("Spatial validation: ",d.blocks.status)
    println("Joint likelihood remains blocked; no production or class output changed.")
    return d
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && acquisition_main()
